#ifndef KKTECQFUNC_H
#define KKTECQFUNC_H
#include <string>
#include <iostream>
#include "wyarray.h"
using namespace std;

class KKTEcqFunc {
    public:
        int m, n;
	int nx, nv, ns;
	int ng, na;

        virtual void G_opr (double *, double *) = 0;
        virtual void Forward(double *y, double *x) = 0;
        virtual void Adjoint(double *y, double *x) = 0;

        virtual void Opr(double *y, double *x, int) = 0;
        virtual void KKTfunc(double *, double *, double *x, double *v) = 0;

	virtual void PrintVector(double *x, int m, char *) = 0;
        virtual void PrintMatrix(double *x, int m, int n, char *) = 0;
	virtual void Free () = 0;
};


#endif
