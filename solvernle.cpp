#include <fstream>
#include "solvernle.h"




/* solve nonlinear programming with augmmented lagrange method */
/* -------------------------------------------------
 * min f(x)  s.t. ve(x) = 0
 *                vn(x) >= 0
 *
 * augmented lagrange formulation
 * L(x,s,ce,cn,sigma) = f(x) + ce' * ve + cn' * (vn - s) 
 *      + 0.5 * sigma * [ve' * ve + (vn - s)' * (vn - s)]
 *
 * where s = max(0, vn - cn / sigma)
 * ------------------------------------------------ */



/* -------------------------------------------------
 * algorithm
 *
 * ---------------------------------------------- */



SolverNle::SolverNle (ProAssign *ass, string fo, string log) {
    this->niter = 200;
    this->outer = 100;
    this->inner = 10;

    this->stop = 0;
    this->flag = 0;

    this->lambda = 0.01;
    this->sigma = 1.000;        // sigma = alpha * sigma
    this->alpha = 2.0;
    this->step = 1.0;

    this->epsilon = 1e-06;
    //this->beta = 0.25;
    this->beta = 0.5;

    this->iprint = 0;
    this->ass = ass;
    this->outputFilename = fo;
    this->logFilename = log;

    this->nx = ass->nx;
    this->ne = ass->ne;
    this->ni = ass->ni;

}

double SolverNle::GetInitStep (double *dx, double *dx_old) {
    double s1 = 0.0;
    double s2 = 0.0;
    for (int j = 0; j < nx; j ++) {
	s1 += dx_old[j] * dx_old[j];
	s2 += dx[j] * dx[j];
    }
    return 0.01 * s1 / (s2 + 1e-30);
}

void SolverNle::UpdateX(double *x, double *x0, double *dx) {
    for (int i = 0; i < this->nx; i ++) {
	x[i] = x0[i] - this->step * dx[i];
	if (x[i] < 1e-30) x[i] = 1e-30;
    }    
}

void SolverNle::UpdateV(double *ce, double *cn, double *ve, double *vn) {
    for (int i = 0; i < this->ne; i ++)
	ce[i] -= this->sigma * ve[i];

    for (int i = 0; i < this->ni; i ++) {
	cn[i] -= this->sigma * vn[i];
	if (cn[i] < 0.0) cn[i] = 0.0;
    }
}

void SolverNle::AugLag (double *x, double *ce, double *cn) {

    double *x0 = new double[nx];
    double *ve = new double[ne];
    double *vn = new double[ni];
    double *vs = new double[ni];
    double *s = new double[ni];

    double ce_norm, ce_norm_old;
    double cn_norm, cn_norm_old;
    double loss;

    // create an output file
    ofstream fp(logFilename);

    // get initial value 
    ass->GetInitPoint(x, ce, cn);
    memcpy(x0, x, nx * sizeof(double));
    //ass->PrintVector(x, nx, "Init x");

    ass->ComputeConstraint(ve, vn, x);
    ass->ComputeSlash(s, cn, vn, this->sigma);
    VectorSub(vs, vn, s, ni);
    ce_norm_old = norm2(ve, ne);
    cn_norm_old = norm2(vs, ni);

    loss = ass->Loss(x);

    for (int iter = 0; iter < this->niter; iter ++) {

	cout << endl << " ### iter = " << iter << "/" << niter << ", sigma = " << this->sigma << ", loss =" << loss << endl;

        this->GradientDescent (x, ce, cn, ve, vn, s, fp);
        loss = ass->Loss(x);

	if (this->stop == 1) break;

	VectorSub(vs, vn, s, ni);
	ce_norm = norm2(ve, ne);
	cn_norm = norm2(vs, ni);

	cout << " ce_norm = " << ce_norm << ", cn_norm = " << cn_norm << endl;
	if (ce_norm < this->epsilon && cn_norm < this->epsilon) break;

	if (ce_norm > this->beta * ce_norm_old &&
	    cn_norm > this->beta * cn_norm_old)
           this->sigma = this->alpha * this->sigma;

	this->UpdateV(ce, cn, ve, vn);
	this->flag += 1;

	ce_norm_old = ce_norm;
	cn_norm_old = cn_norm;

	ass->printResult(outputFilename, x0, x, iter);


    }

    //ass->PrintVector(x, nx, " Update X");

    fp.close();

    cout << endl;
    cout << " --- final x ---- " << endl;
    
    delete [] x0;
    delete [] ve;
    delete [] vn;
    delete [] vs;
    delete [] s;
}



void SolverNle::GradientDescent (double *x, double *ce, double *cn, double *ve, double *vn, double *s, ofstream &fp) {

    double *dx_old = new double[nx];
    double *x0 = new double[nx];

    double *dx = new double[nx];
    double *dce = new double[ne];
    double *dcn = new double[ni];

    double mis, loss;

    for (int iter = 0; iter < this->outer; iter ++) {

	fp << endl << "   * flag = " << this->flag << "/" << this->niter << ", iter = " << iter << "/" << this->outer << ", sigma = " << this->sigma << endl;

	memcpy(x0, x, nx * sizeof(double));

	ass->ComputeConstraint(ve, vn, x);
	ass->ComputeSlash(s, cn, vn, this->sigma);
	ass->ComputeGradient(dx, x, ce, cn, ve, vn, s, this->sigma);
        mis = ass->LossAug(x, ce, cn, ve, vn, s, this->sigma);
	/*
	cout << " mis=" << mis << endl;
	    for (int i = 0; i < nx; i ++)
		    cout << "i=" << i << ", dx=" << dx[i] << endl;
		    */

	//if (this->flag > 0) {
	//if (this->flag == 2 && iter == 0) {
	if (iprint == 1) {
	   cout << endl << " iter = " << iter << "/" << this->outer << endl;
           ass->PrintVector(dx, 10, "update dx");
           ass->PrintVector(x, 10, "update x (outer)");
           ass->PrintVector(ce, 10, "update ce (outer)");
           ass->PrintVector(cn, ni, "update cn (outer)");
           ass->PrintVector(ve, 10, "update ve (outer)");
           ass->PrintVector(vn, ni, "update vn (outer)");
           ass->PrintVector(s, ni, "update s (outer)");
	}

	if (iter > 0) {
           //this->step = this->GetInitStep(dx, dx_old);
	}
	else
           this->step = this->GetInitStep(dx, x);
	if (this->step > 1.0) this->step = 1.0;

	fp << " * iter = " << iter << "/" << this->outer << ", jter = 0/" << this->inner << ", step = " << this->step << ", aug loss = " << mis << endl;

	for (int jter = 0; jter < this->inner; jter ++) {
            this->UpdateX (x, x0, dx);
	    ass->ComputeConstraint(ve, vn, x);
	    ass->ComputeSlash(s, cn, vn, this->sigma);
            loss = ass->LossAug(x, ce, cn, ve, vn, s, this->sigma); 
	    if (iprint == 1) {
	       cout << endl << " jter = " << jter << "/" << this->inner << endl;
               ass->PrintVector(x, 10, "update x (inner)");
               ass->PrintVector(ce, 10, "update ce (inner)");
               ass->PrintVector(cn, ni, "update cn (inner)");
               ass->PrintVector(ve, 10, "update ve (inner)");
               ass->PrintVector(vn, ni, "update vn (inner)");
               ass->PrintVector(s, ni, "update s (inner)");
	    }
	    fp << " * iter = " << iter << "/" << this->outer << ", jter = " << jter + 1 << "/" << this->inner << ", step = " << this->step << ", aug loss = " << loss << endl;
	    if (loss < mis) {
	       this->step = 2.0 * this->step;
     	       break; 
	    }
	    else 
	       this->step = 0.5 * this->step;
	}

        // ass->PrintVector(x, nx, "update x (outer)");
	if (loss >= mis) {
	   memcpy(x, x0, nx * sizeof(double));
	   //this->stop = 1;
	   //break;
	}
	if (fabs(loss / mis - 1.0) < 1e-09) break;

	memcpy(dx_old, dx, nx * sizeof(double));

    }

	//if (this->flag > 0) {
	if (iprint == 1) {
           ass->PrintVector(x, 10, "update x (outer)");
           ass->PrintVector(ce, 10, "update ce (outer)");
           ass->PrintVector(cn, ni, "update cn (outer)");
           ass->PrintVector(ve, 10, "update ve (outer)");
           ass->PrintVector(vn, ni, "update vn (outer)");
           ass->PrintVector(s, ni, "update s (outer)");
	}

    delete [] dx_old;
    delete [] dx;
    delete [] x0;
    delete [] dce;
    delete [] dcn;
}


