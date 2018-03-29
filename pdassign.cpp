#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <set>
#include <map>
#include "wyarray.h"
#include "pdassign.h"
using namespace std;


/* -------------------------------------------------------------------
 * solve Gx = g
 *     G'Gx = G'g
 *       Ax = b
 * where A=G'G is positive definite matrix
 * -----------------------------------------------------------------*/
PDAssign::PDAssign (int nusr, int nmsg) {

   this->nusr = nusr;
   this->nmsg = nmsg;
   this->nx = nusr * nmsg;
   this->ns = nusr + 2 * nmsg + nx;
   this->nv = this->ns;
   this->nm = this->nv + this->nx + this->ns;

   this->b = new double[nv];
   this->v0 = new double[nv];
   this->pf = new double[nv];
   this->pa = new double[nv];
}

void PDAssign::Init(double *v0, double *b) {
   memcpy(this->b, b, nv * sizeof(double));
   memcpy(this->v0, v0, nv * sizeof(double));
}

  
void PDAssign::Free() {
   delete [] this->b;
   delete [] this->v0;
   delete [] this->pf;
   delete [] this->pa;
}



/* y = A'x */
void PDAssign::Adjoint (double *y, double *x) {
    int n1 = nusr + nmsg;
    int n2 = nusr + nmsg * 2;
    for (int i = 0; i < this->nv; i ++)
	pf[i] = x[i] * v0[i];

    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[k + j] += pf[i] - pf[nusr + j] + pf[n1 + j] - pf[n2 + k + j];
	}
    }
}

/* y = Ax */
void PDAssign::Forward (double *y, double *x) {
    int n1 = nusr + nmsg;
    int n2 = nusr + nmsg * 2;
    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[i] += x[k + j];
	    y[nusr + j] -= x[k + j];
	    y[n1 + j] += x[k + j];
	    y[n2 + k + j] -= x[k + j];
	}
    }
    for (int i = 0; i < nv; i ++)
        y[i] *= v0[i];
}

void PDAssign::op(double *y, double *x) {
    memset(pa, 0, nv * sizeof(double));
    this->Forward(pa, x);
    this->Adjoint(y, pa);
}

void PDAssign::printVector(double *x, int m) {
    for (int i = 0; i < m; i ++)
	cout << " i = " << i << ", " << x[i] << endl;
}
