#ifndef ARRAY_H
#define ARRAY_H

#include <stdio.h>

extern int dmax (int, int);
extern int dmin (int, int);
extern float amax (float a, float b);
extern float amin (float a, float b);


// vector
//extern void vector_init_int   (  int *,   int, int);
//extern void vector_init_float (float *, float, int);
extern   int *vector_init_int   (  int, int);
extern float *vector_init_float (float, int);
extern float *matrix2d_to_vector (float **, int, int);
extern float vector_max (float *, int);
extern float vector_min (float *, int);
extern float vector_amax (float *, int);
extern   int vector_amax_index (float *, int);
extern float vector_sum  (float *, int);
extern float vector_asum (float *, int);
extern float vector_ssum (float *, int);
extern float dot_product (float *, float *, int);

extern float *vector_inv (float *, int);
extern float *vector_scale (float *, float, int);
extern float *vector_add (float *, float *, int);
extern float *vector_sub (float *, float *, int);
extern float *vector_mul (float *, float *, int);

void VectorMaxval (float *, int *, float *, int);
void VectorMinval (float *, int *, float *, int);
void VectorMedval (float *, int *, float *, int);

void VectorSet (float *, float, int);
void VectorAbs (float *, float *, int);
void VectorInv (float *, float *, int);
void VectorSqr (float *, float *, int);
void VectorAdd (float *, float *, float *, int);
void VectorSub (float *, float *, float *, int);
void VectorMul (float *, float *, float *, int);
void VectorMin (float *, float *, float *, int);
void VectorMax (float *, float *, float *, int);
void VectorScale (float *, float *, float, int);
void VectorAdd1 (float *, float *, float *, float, float, int);


// two dimensional matrix
extern void free2di   (  int **, int);
extern void free2df   (float **, int);

extern void matrix2dcpy (float **, float **, int, int);

  int   **matrix2d_init_int (  int, int, int);
float **matrix2d_init_float (float, int, int);
float **matrix2d_init_floatl (float, long int, int);

float **vector_to_matrix2d (float *, int, int);
float **matrix2d_product (float **, float **, int, int);
float matrix2d_max  (float **, int, int);
float matrix2d_amax (float **, int, int);

void Matrix2dMul (float **, float **, float **, int, int, int);
void M2dVecMul (float *, float **, float *, int, int);
void VecM2dMul (float *, float *, float **, int, int);

void MatAdd (float **, float **, float **, int, int);
void MatSub (float **, float **, float **, int, int);
void MatMul (float **, float **, float **, int, int);
void MatMin (float **, float **, float **, int, int);
void MatMax (float **, float **, float **, int, int);
void MatInv (float **, float **, int, int);
void MatAbs (float **, float **, int, int);
void MatSet (float **, float, int, int);
void MatScale (float **, float **, float, int, int);
void MatCopy (float **, float **, int, int);
void MatNorm (float **, int, int);
float MatDot (float **, float **, int, int);
float MatAbsMean (float **, int, int);

void MatAdd1 (float **, float **, float **, float, float, int, int);

float dot2_product(float **x, float **y, int m, int n);

// Three dimensional matrix
float  ***matrix3d_init_float (float, int, int, int);
  int  ***matrix3d_init_int   (  int, int, int, int);
float ****matrix4d_init_float (float, int, int, int, int);

extern void malloc3df (float ***, int, int, int);
extern void free3df   (float ***, int, int);
extern void free3di   (  int ***, int, int);
extern void free4df   (float ****, int, int, int);

extern void Mat3dCopy (float ***, float ***, int, int, int, int, int, int);

#endif
