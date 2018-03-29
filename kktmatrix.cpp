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
#include "kktmatrix.h"
using namespace std;


/* -------------------------------------------------------------------
 * min (1/2)x'Gx+x'c
 * s.t. Ax=b
 *
 * KKT matrix:
 * | G  A' |   |-x|   | c|
 * |       | * |  | = |  |
 * | A  0  |   | v|   |-b|
 *
 * where G may be indefinite
 *
 * initial x: Ax-b=0
 *         r: r=Gx+c
 *
 * AA'v = Ar
 * 
 * size(A) = (m, n)
 * size(G) = (n, n)
 * size(x) = n
 * size(v) = m
   -----------------------------------------------------------------*/
KKTMatrix::KKTMatrix (double *g, double *a, int m, int n) {

   this->m = m;
   this->n = n;
   this->ng = n * n;
   this->na = m * n;
   this->nl = m + n;

   this->g = new double[ng];
   memcpy(this->g, g, ng * sizeof(double));

   this->a = new double[na];
   memcpy(this->a, a, na * sizeof(double));

}

  
void KKTMatrix::Free() {
   delete [] this->g;
   delete [] this->a;
}



/* y = Ax */
void KKTMatrix::Forward (double *y, double *x) {
    for (int i = 0; i < m; i ++) {
	int k = i * n;
	for (int j = 0; j < n; j ++) {
	    y[i] += a[k + j] * x[j];
	}
    }
}

/* y = A'x */
void KKTMatrix::Adjoint (double *y, double *x) {
    for (int i = 0; i < m; i ++) {
	int k = i * n;
	for (int j = 0; j < n; j ++) {
	    y[j] += a[k + j] * x[i];
	}
    }
}

/* ifwd = 1, y = A'Ax */
/* ifwd = 0, y = AA'x */
void KKTMatrix::Opr(double *y, double *x, int flag) {
    if (flag == 1) {
       double *p = new double[m];
       memset(p, 0, m * sizeof(double));
       memset(y, 0, n * sizeof(double));
       this->Forward(p, x);
       this->Adjoint(y, p);
       delete [] p;
    }
    else if (flag == 0) {
       double *p = new double[n];
       memset(p, 0, n * sizeof(double));
       memset(y, 0, m * sizeof(double));
       this->Adjoint(p, x);
       this->Forward(y, p);
       delete [] p;
    }
}


/* y = Gx */
void KKTMatrix::G_opr(double *y, double *x) {
    memset(y, 0, n * sizeof(double));
    for (int i = 0; i < n; i ++) {
	int k = i * n;
	for (int j = 0; j < n; j ++) {
            y[i] += g[k + j] * x[j];
	}
    }
}

/* yc = A'v - Gx 
 * yb = Ax        */
void KKTMatrix::KKTfunc(double *yc, double *yb, double *x, double *v) {
    double *p = new double[n];
    memset(p, 0, n * sizeof(double));
    memset(yb, 0, m * sizeof(double));
    this->Forward(yb, x);
    this->Adjoint(p, v);
    this->G_opr(yc, x);
    VectorSub(yc, p, yc, n);
    delete [] p;
}

void KKTMatrix::PrintVector(double *x, int m, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) 
       cout << " i = " << i << ", " << x[i] << endl;
}
  
void KKTMatrix::PrintMatrix(double *x, int m, int n, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) {
       cout << " i=" << i << " "; 
       for (int j = 0; j < n; j ++)
           cout << x[i*n+j] << " ";
       cout << endl;
   }
}
  
