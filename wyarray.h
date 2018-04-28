#ifndef WYARRAY_H
#define WYARRAY_H

#include <stdio.h>



double *vector_init_double (double a, int n);
double **matrix2d_init_double (double a, int m, int n);

double norm2 (double *x, int n);


int vector_max_index (double *x, int n);
double vector_min (double *x, int n);
double vector_max (double *x, int n);
double vector_amax (double *, int);
double vector_sum (double *a, int n);
double dot_product (double *, double *, int);
double dot2_product(double **x, double **y, int m, int n);

void VectorAdd (double *z, double *x, double *y, double a, double b, int n);
void VectorSub (double *y, double *a, double *b, int n);
void M2dVecMul (double *y, double **a, double *b, int m, int n);
#endif
