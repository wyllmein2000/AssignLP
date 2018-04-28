#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "inout.h"
#include "wyarray.h"
#include "proassign.h"
using namespace std;


class SolverAssDual {
    private:
        int niter;
        int inner;
	int iprint;
	int nx, nusr, nmsg;

	double epsilon;
	double *w, *b;

	string logFilename;
	string outFilename;

    public:
        SolverAssDual (double *, double *, int, int);
	void SetFilename(string, string);

	void GetLambdaFromMu (int *, double *, double *);
	void GetInitDualValue (double *, double *, int *);
	void ComputeGradient (double *);

	double GetInitStep (double *, double *);
	double DualLoss(double *, double *);
	double PrimalLoss(double *);

	void UpdateMu (double *, double *, double *, double);
	void UpdateX(double *, int *);
	void GradientDescent(double *);

};
