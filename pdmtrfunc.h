#ifndef PDMTRFUNC_H
#define PDMTRFUNC_H
#include <string>
#include <iostream>
using namespace std;

class PDMtrFunc {
    private:
	double *g;
    public:
	int nx, ny, nm;

	virtual void Free () = 0;
	virtual void Init (double *, int, int) = 0;
	virtual void Forward (double *, double *) = 0;
	virtual void Adjoint (double *, double *) = 0;
	virtual void Opr (double *, double *) = 0;
};


#endif
