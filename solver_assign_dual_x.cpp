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
    this->niter = 200;
    this->inner = 10;
    this->epsilon = 1e-09;

    this->nusr = nusr;
    this->nmsg = nmsg;
    this->nx = nusr * nmsg;
    this->w = new double[this->nx];
    this->b = new double[this->nmsg];

    memcpy(this->w, w, nx * sizeof(double));
    memcpy(this->b, b, nmsg * sizeof(double));

    this->iprint = 0;
    this->logFilename = "output/log.txt";
    this->outFilename = "output/out.txt";
}

void SolverAssDual::SetFilename(string fo, string log) {
    this->outFilename = fo;
    this->logFilename = log;
}

void SolverAssDual::GetLambdaFromMu(int *index, double *y, double *x) {
    for (int i = 0; i < this->nusr; i ++) {
	int k = i * this->nmsg;
	double a = 0.0;
	y[i] = 0.0;
	index[i] = 0;
	for (int j = 0; j < this->nmsg; j ++) {
	    a = this->w[k + j] + x[j];
            if (y[i] < a) {
	       y[i] = a;
	       index[i] = j;
	    }
	}
    }
}


void SolverAssDual::GetInitDualValue(double *lambda, double *mu, int *index) {
    memset(mu, 0, nmsg * sizeof(double));
    this->GetLambdaFromMu(index, lambda, mu);
}

void SolverAssDual::ComputeGradient(double *grad) {
    for (int j = 0; j < nmsg; j ++)
	grad[j] = -1.0 * this->b[j];
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

double SolverAssDual::DualLoss(double *lambda, double *mu) {
    double y = 0.0;
    for (int i = 0; i < this->nusr; i ++)
	y += lambda[i];
    for (int j = 0; j < this->nmsg; j ++)
	y -= mu[j] * this->b[j];
    return y;
}

double SolverAssDual::PrimalLoss(double *x) { 
    int k;
    double y = 0.0;
    for (int i = 0; i < nusr; i ++) {
        k = i * nmsg;
        for (int j = 0; j < nmsg; j ++)
	    y += x[k + j] * this->w[k + j];
    }
    return y;
}

void SolverAssDual::UpdateMu(double *mu, double *mu0, double *grad, double step) {
    for (int j = 0; j < nmsg; j ++) {
	mu[j] = mu0[j] - step * grad[j];
	if (mu[j] < 0.0) mu[j] = 0.0;
    }
}

void SolverAssDual::UpdateX(double *x, int *index) {
    int j;
    memset(x, 0, nusr * nmsg * sizeof(double));
    for (int i = 0; i < nusr; i ++) {
	j = index[i];
	x[i * nmsg + j] = 1.0;
    }
}

void SolverAssDual::GradientDescent(double *x) {
    
    double *x0 = new double[nx];
    double *lambda = new double[nusr];
    double *mu0 = new double[nmsg];
    double *mu = new double[nmsg];
	    
    double *grad = new double[nmsg];
    double *gold = new double[nmsg];
    int *jmax = new int[nusr];

    double loss, mis, mis_pri;
    double step = 1.0;
    double loss_perc = 1.0;

    ProAssign *ass = new ProAssign(this->w, this->b, nusr, nmsg);


    // create an output file
    ofstream fp(logFilename);

    // set initial value
    this->GetInitDualValue(lambda, mu, jmax);
    this->UpdateX(x0, jmax);

    // initial loss
    mis = this->DualLoss(lambda, mu);
    mis_pri = this->PrimalLoss(x0);

    // ???
    this->ComputeGradient(grad);

    for (int iter = 0; iter < this->niter; iter ++) {

	cout << endl << " ### iter = " << iter << "/" << niter << ", perc = " << loss_perc << ", primal loss =" << mis_pri << ", dual loss =" << mis << endl;
	fp << endl << " ### iter = " << iter << "/" << niter << ", perc = " << loss_perc << ", primal loss =" << mis_pri << ", dual loss =" << mis << endl;

	memcpy(gold, grad, nmsg * sizeof(double));
	memcpy(mu0, mu, nmsg * sizeof(double));

	this->ComputeGradient(grad);
	for (int i = 0; i < nmsg; i ++)
	    cout << " i = " << gold[i] << ", " << grad[i] << endl;

	if (iter == 0) {
	    step = this->GetInitStep(grad, gold);
	    if (step > 1.0) step = 1.0;
	}
	step = 1e-05;
	cout << " step = " << step << endl;

	for (int jter = 0; jter < this->inner; jter ++) {

            this->UpdateMu (mu, mu0, grad, step);
	    this->GetLambdaFromMu (jmax, lambda, mu);
            loss = this->DualLoss (lambda, mu);

	    cout << " * iter = " << iter << "/" << this->niter << ", jter = " << jter + 1 << "/" << this->inner << ", step = " << step << ", dual loss = " << loss << endl;
	    fp << " * iter = " << iter << "/" << this->niter << ", jter = " << jter + 1 << "/" << this->inner << ", step = " << step << ", dual loss = " << loss << endl;

	    if (loss < mis) {
	       step = 2.0 * step;
     	       break; 
	    }
	    else 
	       step = 0.5 * step;
	}

	//if (loss >= mis) break;

	this->UpdateX(x, jmax);
        mis = this->DualLoss(lambda, mu);
	mis_pri = this->PrimalLoss(x);
	ass->printResult(outFilename, x0, x, iter);
	if (loss >= mis) break;

	loss_perc = fabs(loss / mis - 1.0);
	if (loss_perc < this->epsilon) break;


    }

    delete [] gold;
    delete [] grad;
    delete [] jmax;
    delete [] x0;
    delete [] mu;
    delete [] mu0;
    delete [] lambda;
}


