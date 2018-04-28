#include <fstream>
#include "solver_assign_dual.h"




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



SolverAssDual::SolverAssDual (double *w, double *b, int nusr, int nmsg) {
    this->niter = 500;
    this->inner = 10;
    this->epsilon = 1e-09;
    this->lambda = 0.001;

    this->nusr = nusr;
    this->nmsg = nmsg;
    this->nx = nusr * nmsg;
    this->w = new double[this->nx];
    this->b = new double[this->nmsg];

    memcpy(this->w, w, nx * sizeof(double));
    memcpy(this->b, b, nmsg * sizeof(double));

    this->cm = new double[this->nusr];

    this->iprint = 0;
    this->logFilename = "output/log.txt";
    this->outFilename = "output/out.txt";
}

void SolverAssDual::SetFilename(string fo, string log) {
    this->outFilename = fo;
    this->logFilename = log;
}

void SolverAssDual::GetInitDualValue(double *alpha) {
    memset(alpha, 0, this->nmsg * sizeof(double));
}

void SolverAssDual::UpdatePrimal(double *x, double *z, double *alpha) {
    for (int i = 0; i < this->nusr; i ++) {
	int k = i * this->nmsg;
	z[i] = 0.0;

	double cc = 0.0;
        for (int j = 0; j < this->nmsg; j ++) {
	    double c1 = (this->w[k + j] + alpha[j]) / this->lambda;
	    if (cc < c1) cc = c1;
	}
	this->cm[i] = cc;

        for (int j = 0; j < this->nmsg; j ++) {
	    x[k + j] = exp((this->w[k + j] + alpha[j]) / this->lambda - cc);
	    z[i] += x[k + j];
	}
        for (int j = 0; j < this->nmsg; j ++)
	    x[k + j] /= z[i];
    }
}

void SolverAssDual::ComputeGradient(double *grad, double *x) {
    for (int j = 0; j < nmsg; j ++) {
	grad[j] = -1.0 * this->b[j];
	for (int i = 0; i < nusr; i ++)
	    grad[j] += x[i * nmsg + j];
    }
}

double SolverAssDual::GetInitStep(double *grad, double *gold) {
    double s1 = 0.0;
    double s2 = 0.0;
    for (int j = 0; j < nmsg; j ++) {
	s1 += gold[j] * gold[j];
	s2 += grad[j] * grad[j];
    }
    return 0.01 * s1 / (s2 + 1e-30);
}

double SolverAssDual::LossDual(double *z, double *alpha) {
    double y = 0.0;
    for (int i = 0; i < this->nusr; i ++)
	y += log(z[i]) + this->cm[i];
    y *= this->lambda;
    for (int j = 0; j < this->nmsg; j ++)
	y -= alpha[j] * this->b[j];
    return y;
}

double SolverAssDual::LossPrim(double *x) { 
    int k, m;
    double y = 0.0;
    for (int i = 0; i < nusr; i ++) {
        k = i * nmsg;
        for (int j = 0; j < nmsg; j ++) {
            m = k + j;
	    if (x[m] > 1e-30)
	    y += x[m] * (this->w[m] - this->lambda * log(x[m]));
	    //cout << i << "," << j << " " << x[m] << " " << log(x[m]) << endl;
	}
    }
    return y;
}

void SolverAssDual::UpdateAlpha(double *mu, double *mu0, double *grad, double step) {
    for (int j = 0; j < nmsg; j ++) {
	mu[j] = mu0[j] - step * grad[j];
	if (mu[j] < 0.0) mu[j] = 0.0;
    }
}

void SolverAssDual::GradientDescent(double *x) {
    
    double *alpha0 = new double[nmsg];
    double *alpha = new double[nmsg];
    double *x0 = new double[nx];
    double *z = new double[nusr];

    double *grad = new double[nmsg];
    double *gold = new double[nmsg];

    double loss, misDual, misPrim;
    double step = 1.0;
    double loss_perc = 1.0;

    ProAssign *ass = new ProAssign(this->lambda, this->w, this->b, nusr, nmsg, "entropy");

    //for (int i = 0; i < nmsg; i ++)
//	    cout << "i=" << i << ", " << this->b[i] << endl;

    // create an output file
    ofstream fp(logFilename);

    // set initial value
    this->GetInitDualValue(alpha);
    this->UpdatePrimal(x0, z, alpha);
    memcpy(x, x0, nx * sizeof(double));

    // initial loss
    misDual = this->LossDual(z, alpha);
    misPrim = this->LossPrim(x0);

    // ???
    this->ComputeGradient(grad, x0);

    for (int iter = 0; iter < this->niter; iter ++) {

	if (iter % 10 == 0)
	   iprint = 1;
	else
	   iprint = 0;

	if (iter % 10 == 0) {
	fp << endl << " ### iter = " << iter << "/" << niter << ", perc = " << loss_perc << ", primal loss =" << misPrim << ", dual loss =" << misDual << endl;
	}
	if (iprint == 1)
	cout << endl << " ### iter = " << iter << "/" << niter << ", perc = " << loss_perc << ", primal loss =" << misPrim << ", dual loss =" << misDual << endl;

	memcpy(alpha0, alpha, nmsg * sizeof(double));
	memcpy(gold, grad, nmsg * sizeof(double));

	this->ComputeGradient(grad, x);
	//for (int i = 0; i < nmsg; i ++)
	//    cout << " i = " << gold[i] << ", " << grad[i] << endl;

	if (iter == 0) {
	    //step = this->GetInitStep(grad, gold);
	    step = 0.1 / nmsg;
	    if (step > 1.0) step = 1.0;
	}

	for (int jter = 0; jter < this->inner; jter ++) {

            this->UpdateAlpha (alpha, alpha0, grad, step);
	    this->UpdatePrimal (x, z, alpha);
            loss = this->LossDual (z, alpha);
	    /*
	    for (int i = 0; i < nmsg; i ++)
		    cout << i << " alpha,alpha0,grad=" << alpha[i] << " " << alpha0[i] << " " << grad[i] << endl;
	    for (int i = 0; i < nusr; i ++)
		    cout << i << " z=" << z[i] << endl;
	    exit(0); */

	    if (iter % 10 == 0) {
	    fp << " * iter = " << iter << "/" << this->niter << ", jter = " << jter + 1 << "/" << this->inner << ", step = " << step << ", dual loss = " << loss << endl;
	    }
	    if (iprint == 1)
	    cout << " * iter = " << iter << "/" << this->niter << ", jter = " << jter + 1 << "/" << this->inner << ", step = " << step << ", dual loss = " << loss << endl;

	    if (loss < misDual) {
	       step = 2.0 * step;
     	       break; 
	    }
	    else 
	       step = 0.5 * step;
	}

	if (iter % 100 == 0) 
	   ass->printResult(outFilename, x0, x, iter);

	loss_perc = fabs(loss / misDual - 1.0);
	if (loss >= misDual || loss_perc < this->epsilon) {
           ass->printResult(outFilename, x0, x, iter);
	   break;
	}

        misPrim = this->LossPrim(x);
        misDual = loss;
    }


    delete [] gold;
    delete [] grad;
    delete [] x0;
    delete [] z;
    delete [] alpha;
    delete [] alpha0;
}


