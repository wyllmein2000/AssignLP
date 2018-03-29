#ifndef PDASSIGN_H
#define PDASSIGN_H
#include <string>
#include <iostream>
using namespace std;

class PDAssign: public PDMtrFunc {
    private:
	int nusr;
	int nmsg;
	int ns, nx, nv;
	int nm;
	double *v0, *b;
	double *pf, *pa;

    public:

        PDAssign (int nusr, int nmsg);

	void Free ();
	void Init (double *, double *);
        void Forward(double *y, double *x);
        void Adjoint(double *y, double *x);
        void op(double *y, double *x);

        void printVector(double *x, int m);
};


#endif
