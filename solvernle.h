#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "inout.h"
#include "wyarray.h"
#include "proassign.h"
using namespace std;


class SolverNle {
    private:
        int niter;
        int inner;
	int outer;
	int iprint;
	int flag;
	int stop;

	double lambda;
	double sigma, alpha;
	double epsilon, beta;
	double step;

	string logFilename;
	string outputFilename;

    public:
	int nx, ne, ni;
	ProAssign *ass;

        SolverNle (ProAssign *, string, string);
	double GetInitStep(double *, double *);

	void UpdateX(double *, double *, double *);
	void UpdateV(double *, double *, double *, double *);
	void AugLag (double *, double *, double *);
	void GradientDescent(double *, double *, double *, double *, double *, double *, ofstream &);

};
