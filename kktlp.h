#ifndef KKTLP_H
#define KKTLP_H
#include <string>
#include <iostream>
//#include "kktecqfunc.h"
#include "kktmatrix.h"
using namespace std;

//class KKTlp: public KKTEcqFunc {
//class KKTlp: public KKTMatrix {
class KKTlp {

    public:
	int m, n, na, ng;
        double *a, *g;              // const parameters
	double *x0, *v0, *s0;       // parameters may change in iteration

        KKTlp (double *, int, int);

	// need to self-define
	void G_opr (double *, double *);
        void Forward(double *, double *);
        void Adjoint(double *, double *);
	void Precon(double *, double *, double *, int);

        void Opr (double *, double *, int);
        void KKTfunc(double *, double *, double *, double *);

	void PrintVector(double *x, int m, char *);
        void PrintMatrix(double *x, int m, int n, char *str);
	void UpdateX (double *, double *, double *);
	void Free ();
};
#endif
