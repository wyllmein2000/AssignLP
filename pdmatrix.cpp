#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "pdmatrix.h"
#include "wyarray.h"
#include "inout.h"
using namespace std;

/*
 *   G x = g
 * G'G x = G'g
 *   A x = b
 *
 * where G is ny x nx matrix
 */

PDMatrix::PDMatrix (double *g, int nx, int ny) {
    this->nx = nx;
    this->ny = ny;
    this->nm = ny;
    this->n = nx * ny;
    this->g = new double[n];
    memcpy(this->g, g, n * sizeof(double));
}

void PDMatrix::Init (double *x, int nx, int ny) {
    this->nx = nx;
    this->ny = ny;
    this->nm = ny;
    this->n = nx * ny;

    this->g = new double[n];
    memcpy(this->g, x, n * sizeof(double));
}

void PDMatrix::Forward(double *y, double *x) {
   int k;
   for (int i=0; i<nx; i++) {
       k=i*nx;
       for (int j=0; j<ny; j++)
           y[i] += this->g[k+j]*x[j];
   }
}

void PDMatrix::Adjoint(double *y, double *x) {
   int k;
   for (int i=0; i<nx; i++) {
       k=i*nx;
       for (int j=0; j<ny; j++)
           y[j] += this->g[k+j]*x[i];
   }
}

void PDMatrix::Opr (double *y, double *x) {
   double *p = new double[nx];
   memset(p, 0, nx * sizeof(double));
   memset(y, 0, ny * sizeof(double));
   this->Forward(p, x);
   this->Adjoint(y, p);
   delete [] p;
}

void PDMatrix::Free () {
   delete [] this->g;
}


void PDMatrix::PrintVector(double *x, int m, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) 
       cout << " i = " << i << ", " << x[i] << endl;
}
