#ifndef KKTMATRIX_H
#define KKTMATRIX_H
#include <string>
#include <iostream>
//#include "kktlp.h"
//#include "kktecqfunc.h"
//#include "kktmatrix.h"
using namespace std;

//class KKTMatrix: public KKTlp {
//class KKTMatrix: public KKTEcqFunc {
class KKTMatrix {

    public:
	int m, n, ng, na, nl;
        double *g, *a;

        KKTMatrix (double *, double *, int, int);
        KKTMatrix (double *, int, int);

	void G_opr (double *, double *);
        void Forward(double *y, double *x);
        void Adjoint(double *y, double *x);
        void Opr (double *y, double *x, int);

        void KKTfunc(double *, double *, double *x, double *v);

	void PrintVector(double *x, int m, char *);
        void PrintMatrix(double *x, int m, int n, char *str);
	void Free ();
};
#endif
