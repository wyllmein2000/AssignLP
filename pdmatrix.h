#ifndef PDMatrix_H
#define PDMatrix_H
#include <string>
#include <iostream>
#include "pdmtrfunc.h"
using namespace std;

class PDMatrix: public PDMtrFunc {
    private:
        double *g;
    public:
	int nx, ny, n;
	int nm;

        PDMatrix (double *, int, int);

	void Free ();
	void Init(double *, int, int);

        void Forward(double *y, double *x);
        void Adjoint(double *y, double *x);
        void Opr(double *y, double *x);

	void PrintVector(double *x, int m, char *);
};


#endif
