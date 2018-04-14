#ifndef KKTECQFUNC_H
#define KKTECQFUNC_H
#include <string>
#include <iostream>
#include "wyarray.h"
#include "pdmtrfunc.h"
#include "cglv.h"
using namespace std;

class KKTEcqFunc {
    private:
	int nusr;
	int nmsg;
	int ns, nx, nv;
	int nm;

        double lambda;
	double gamma;

	double *w, *c;
	double *x0, *v0, *s0;
	double *pf, *pa;
	double *xo;

    public:

        KKTEcqFunc (double *score, double *bottom, double *upper, int nusr, int nmsg);

	void init ();
	void free ();
        void GetX (double *x);
        void GetDim (int *);
        void UpdateX (double *x);
	void GetVo ();

        void Forward(double *y, double *x);
        void Adjoint(double *y, double *x);
        void op(double *y, double *x);
        void G_opr (double *, double *);
        void P_opr (double *, double *);

        void residual_aff(double *y, double *x);
        void maxres(double *r1, double *r2, double *x);

        void round(double *y, double *x);
        void flow(double *y, double *x);
        double entropy(double *x, int n);
        double yield(double *x);
        double misfit(double *x);

        void printResult(string outputFileName);
}
;


#endif
