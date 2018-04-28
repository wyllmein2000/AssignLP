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

	double lambda;
	double epsilon;
	double *w, *b;
	double *cm;

	string logFilename;
	string outFilename;

    public:
        SolverAssDual (double *, double *, int, int);
	void SetFilename(string, string);

	void GetInitDualValue (double *);
	void UpdatePrimal(double *, double *, double *);
	void ComputeGradient (double *, double *);

	double GetInitStep (double *, double *);
	double LossDual(double *, double *);
	double LossPrim(double *);

	void UpdateAlpha (double *, double *, double *, double);
	void GradientDescent(double *);

};
