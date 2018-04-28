#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wyarray.h"

    double *vector_init_double (double a, int n) {
      double *x=(double *)malloc(n*sizeof(double));
      for (int i=0; i<n; i++)
          x[i]=a;
      return x;
   }


  double **matrix2d_init_double (double a, int m, int n) {
  double **x=(double **)malloc(m*sizeof(double *));
  // x[0]=(double *)malloc(sizeof(double)*n*m);
  for (int i=0; i<m; i++)
      x[i]=(double *)malloc(sizeof(double)*n);
      // x[i]=x[i]+n;

  for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++)
          x[i][j]=a;
  }
  return x;
  }


   void M2dVecMul (double *y, double **a, double *b, int m, int n) {
   for (int j=0; j<m; j++) {
       y[j]=0.0;
       for (int k=0; k<n; k++)
           y[j] += a[j][k]*b[k];
   }
}


   void VectorAdd (double *z, double *x, double *y, double a, double b, int n) {
      for (int i=0; i<n; i++) z[i]=a*x[i]+b*y[i];
   }


   void VectorSub (double *y, double *a, double *b, int n) {
      for (int i=0; i<n; i++) y[i]=a[i]-b[i];
   }


   double vector_sum (double *a, int n) {
      double y = 0.0;
      for (int i=0; i<n; i++)
          y += a[i];
      return y;
   }

// ===================================================================
   double vector_min (double *x, int n) {
   double y=x[0];
   for (int i=1; i<n; i++)
       if (y>x[i]) y=x[i];
   return y;
}
// ===================================================================
   double vector_max (double *x, int n) {
   double y=x[0];
   for (int i=1; i<n; i++)
       if (y<x[i]) y=x[i];
   return y;
}
// ===================================================================
   double vector_amax (double *x, int n) {
   double y=x[0];
   for (int i=1; i<n; i++)
       if (y<fabs(x[i])) y=fabs(x[i]);
   return y;
}
// ===================================================================
   int vector_max_index (double *x, int n) {
   int k=0;
   double y=x[0];
   for (int i=1; i<n; i++) {
       if (y<x[i]) {
          y=x[i];
          k=i;
       }
   }
   return k;
}
// ===================================================================
   double dot_product (double *a, double *b, int n) {
      double y = 0.0;
      for (int i=0; i<n; i++)
          y += a[i]*b[i];
      return y;
   }
// ===================================================================

// ===================================================================
   double norm2 (double *x, int n) {
      double y = 0.0;
      for (int i=0; i<n; i++)
          y += x[i]*x[i];
      return sqrt(y)/n;
   }
// ===================================================================
