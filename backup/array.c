#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ===================================================================
   int dmax (int a, int b) {
       int c=a>b?a:b;
       return c;
   }
// ===================================================================
   int dmin (int a, int b) {
       int c=a>b?b:a;
       return c;
   }
// ===================================================================
   float amax (float a, float b) {
       float c=a>b?a:b;
       return c;
   }
// ===================================================================
   float amin (float a, float b) {
       float c=a>b?b:a;
       return c;
   }
// ===================================================================





// operation on 1d vector
// ===================================================================
   int *vector_init_int (int a, int n) {
      int *x=(int *)malloc(n*sizeof(int));
      for (int i=0; i<n; i++)
          x[i]=a;
      return x;
   }
// ===================================================================
   float *vector_init_float (float a, int n) {
      float *x=(float *)malloc(n*sizeof(float));
      for (int i=0; i<n; i++)
          x[i]=a;
      return x;
   }
// ===================================================================
   float vector_max (float *x, int n) {
   float y=x[0];
   for (int i=1; i<n; i++)
       if (y<x[i]) y=x[i];
   return y;
}
// ===================================================================
   float vector_min (float *x, int n) {
   float y=x[0];
   for (int i=1; i<n; i++)
       if (y>x[i]) y=x[i];
   return y;
}
// ===================================================================
   float vector_amax (float *x, int n) {
   float y=x[0];
   for (int i=1; i<n; i++)
       if (y<fabs(x[i])) y=fabs(x[i]);
   return y;
}
// ===================================================================
   int vector_amax_index (float *x, int n) {
   int k=0;
   float y=x[0];
   for (int i=1; i<n; i++) {
       if (y<fabs(x[i])) {
          y=fabs(x[i]);
          k=i;
       }
   }
   return k;
}
// ===================================================================
   float vector_sum (float *a, int n) {
      float y = 0.0;
      for (int i=0; i<n; i++)
          y += a[i];
      return y;
   }
// ===================================================================
   float vector_asum (float *a, int n) {
      float y = 0.0;
      for (int i=0; i<n; i++)
          y += fabs(a[i]);
      return y;
   }
// ===================================================================
   float vector_ssum (float *a, int n) {
      float y = 0.0;
      for (int i=0; i<n; i++)
          y += a[i]*a[i];
      return y;
   }
// ===================================================================
   float dot_product (float *a, float *b, int n) {
      float y = 0.0;
      for (int i=0; i<n; i++)
          y += a[i]*b[i];
      return y;
   }
// ===================================================================
   float *vector_inv (float *x, int n) {
   float *y=(float *)malloc(n*sizeof(float));
   for (int i=0; i<n; i++) {
       if (x[i]==0.0)
          y[i]=0.0;
       else
          y[i]=1.0/x[i];
   }
   return y;
}
// ===================================================================
   float *vector_scale (float *x, float a, int n) {
   float *y=(float *)malloc(n*sizeof(float));
   for (int i=0; i<n; i++) y[i]=x[i]*a;
   return y;
}
// ===================================================================
   float *vector_add (float *a, float *b, int n) {
      float *y=(float *)malloc(n*sizeof(float));
      for (int i=0; i<n; i++) y[i]=a[i]+b[i];
      return y;
   }
// ===================================================================
   float *vector_sub (float *a, float *b, int n) {
      float *y=(float *)malloc(n*sizeof(float));
      for (int i=0; i<n; i++) y[i]=a[i]-b[i];
      return y;
   }
// ===================================================================
   float *vector_mul (float *a, float *b, int n) {
      float *y=(float *)malloc(n*sizeof(float));
      for (int i=0; i<n; i++) y[i]=a[i]*b[i];
      return y;
   }
// ===================================================================





// ===================================================================
   void VectorMaxval (float *a, int *k, float *x, int n) {
   int m=0;
   float y=x[0];
   for (int i=1; i<n; i++) {
       if (y<x[i]) {
          y=x[i];
          m=i;
       }
   }
   *a=y;
   *k=m;
   }
// ===================================================================
   void VectorMinval (float *a, int *k, float *x, int n) {
   int m=0;
   float y=x[0];
   for (int i=1; i<n; i++) {
       if (y>x[i]) {
          y=x[i];
          m=i;
       }
   }
   *a=y;
   *k=m;
   }
// ===================================================================
   void VectorMedval (float *a, int *k, float *x, int n) {
   int *m=(int *)malloc(n*sizeof(int));
   float *y=(float *)malloc(n*sizeof(float));

   memcpy(y,x,n*sizeof(float));
   for (int i=0; i<n; i++) m[i]=i;

   int nh=n/2;

   for (int i=0; i<nh; i++) {
       int km=i;
       float xm=y[i];
       for (int j=i+1; j<n; j++) {
           if (xm>y[j]) {
              xm=y[j];
              km=j;
           }
       }
       if (i!=km) {
          y[km]=y[i];
          y[i]=xm;
          m[km]=i;
          m[i]=km;
       }
   }

   *a=y[nh-1];
   *k=m[nh-1];

   free(y);
   free(m);
   }
// ===================================================================
   void VectorSet (float *x, float a, int n) {
      for (int i=0; i<n; i++) x[i]=a;
   }
// ===================================================================
   void VectorAbs (float *y, float *x, int n) {
   for (int i=0; i<n; i++) y[i]=fabs(x[i]);
}
// ===================================================================
   void VectorSqr (float *y, float *x, int n) {
   for (int i=0; i<n; i++) 
       if (x[i]<0)
          y[i]=0.0;
       else
          y[i]=sqrt(x[i]);
}
// ===================================================================
   void VectorInv (float *y, float *x, int n) {
      for (int i=0; i<n; i++) {
          if (x[i]==0.0)
             y[i]=0.0;
          else
             y[i]=1.0/x[i];
      }
   }
// ===================================================================
   void VectorAdd (float *y, float *a, float *b, int n) {
      for (int i=0; i<n; i++) y[i]=a[i]+b[i];
   }
// ===================================================================
   void VectorSub (float *y, float *a, float *b, int n) {
      for (int i=0; i<n; i++) y[i]=a[i]-b[i];
   }
// ===================================================================
   void VectorMul (float *y, float *a, float *b, int n) {
      for (int i=0; i<n; i++) y[i]=a[i]*b[i];
   }
// ===================================================================
   void VectorScale (float *y, float *x, float a, int n) {
      for (int i=0; i<n; i++) y[i]=x[i]*a;
   }
// ===================================================================
   void VectorMin (float *y, float *a, float *b, int n) {
      for (int i=0; i<n; i++) {
          if (a[i]<b[i])
             y[i]=a[i];
          else
             y[i]=b[i];
      }
   }
// ===================================================================
   void VectorMax (float *y, float *a, float *b, int n) {
      for (int i=0; i<n; i++) {
          if (a[i]<b[i])
             y[i]=b[i];
          else
             y[i]=a[i];
      }
   }
// ===================================================================
   void VectorAdd1 (float *z, float *x, float *y, float a, float b, int n) {
      for (int i=0; i<n; i++) z[i]=a*x[i]+b*y[i];
   }
// ===================================================================









// operation on 2d or multi-dimension matrix
// ===================================================================
  int **matrix2d_init_int (int a, int m, int n) {
  int **x=(int **)malloc(m*sizeof(int *));
  for (int i=0; i<m; i++)
      x[i]=(int *)malloc(sizeof(int)*n);

  for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++)
          x[i][j]=a;
  }
  return x;
}
// ===================================================================
  float **matrix2d_init_float (float a, int m, int n) {
  float **x=(float **)malloc(m*sizeof(float *));
  // x[0]=(float *)malloc(sizeof(float)*n*m);
  for (int i=0; i<m; i++)
      x[i]=(float *)malloc(sizeof(float)*n);
      // x[i]=x[i]+n;

  for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++)
          x[i][j]=a;
  }
  return x;
}
// ===================================================================
  float **matrix2d_init_floatl (float a, long int m, int n) {
  float **x=(float **)malloc(m*sizeof(float *));
  for (long int i=0; i<m; i++)
      x[i]=(float *)malloc(sizeof(float)*n);

  for (long int i=0; i<m; i++) {
      for (int j=0; j<n; j++)
          x[i][j]=a;
  }
  return x;
}
// ===================================================================


// ===================================================================
// ===================================================================
   float ***matrix3d_init_float (float a, int nx, int ny, int nz) {
   float ***x=(float ***)malloc(nx*sizeof(float **));

   for (int i=0; i<nx; i++)
       x[i]=(float **)malloc(ny*sizeof(float *));

   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++)
       x[i][j]=(float *)malloc(nz*sizeof(float));
   }
   
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
   for (int k=0; k<nz; k++)
       x[i][j][k]=a;
   }
   }
   return x;
}
// ===================================================================
   int ***matrix3d_init_int (int a, int nx, int ny, int nz) {
   int ***x=(int ***)malloc(nx*sizeof(int **));

   for (int i=0; i<nx; i++)
       x[i]=(int **)malloc(ny*sizeof(int *));

   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++)
       x[i][j]=(int *)malloc(nz*sizeof(int));
   }
   
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
   for (int k=0; k<nz; k++)
       x[i][j][k]=a;
   }
   }
   return x;
}
// ===================================================================
   float ****matrix4d_init_float (float a, int nx, int ny, int nz, int m) {
   float ****x=(float ****)malloc(nx*sizeof(float ***));

   for (int i=0; i<nx; i++)
       x[i]=(float ***)malloc(ny*sizeof(float **));

   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++)
       x[i][j]=(float **)malloc(nz*sizeof(float *));
   }
   
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
   for (int k=0; k<nz; k++)
       x[i][j][k]=(float *)malloc(m*sizeof(float));
   }
   }

   for (int i=0; i<nx; i++)
   for (int j=0; j<ny; j++)
   for (int k=0; k<nz; k++)
   for (int l=0; l< m; l++)
       x[i][j][k][l]=a;

   return x;
}
// ===================================================================




// ===================================================================
void free2di (int **x, int m) {
   // free(x[0]);
   // free(x);
   for (int i=0; i<m; i++)
      free(x[i]);
  free(x);
}
// ===================================================================
void free2df (float **x, int m) {
  // free(x[0]);
  // free(x);
  for (int i=0; i<m; i++)
      free(x[i]);
  free(x);
}
// ===================================================================
void free3df (float ***x, int nx, int ny) {
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
       free(x[i][j]);
       x[i][j]=0;
   }
   free(x[i]);
   x[i]=0;
   }
   free(x);
}
// ===================================================================
void free3di (int ***x, int nx, int ny) {
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
       free(x[i][j]);
       x[i][j]=0;
   }
   free(x[i]);
   x[i]=0;
   }
   free(x);
}
// ===================================================================
   void free4df (float ****x, int nx, int ny, int nz) {
   for (int i=0; i<nx; i++) {
   for (int j=0; j<ny; j++) {
   for (int k=0; k<nz; k++) {
       free(x[i][j][k]);
       x[i][j][k]=0;
   }
       free(x[i][j]);
       x[i][j]=0;
   }
   free(x[i]);
   x[i]=0;
   }
   free(x);
}
// ===================================================================















// ===================================================================
   float dot2_product(float **x, float **y, int m, int n) {
   float z=0.0;
   for (int i=0; i<m; i++)
   for (int j=0; j<n; j++)
       z += x[i][j]*y[i][j];
   return z;
}
// ===================================================================
   void matrix2dcpy (float **y, float **x, int m, int n) {
   for (int j=0; j<m; j++) {
       for (int i=0; i<n; i++) 
           y[j][i]=x[j][i];
   }
}
// ===================================================================
   float **vector_to_matrix2d (float *x, int m, int n) {
   float **y=matrix2d_init_float(0.0,m,n);
   for (int j=0; j<m; j++) {
       int k=j*n;
   for (int i=0; i<n; i++)
       y[j][i]=x[k+i];
   }
   return y;
}
// ===================================================================
   float *matrix2d_to_vector (float **x, int m, int n) {
   float *y=(float *)malloc(m*n*sizeof(float));
   for (int j=0; j<m; j++) {
       int k=j*n;
   for (int i=0; i<n; i++)
       y[k+i]=x[j][i];
   }
   return y;
}
// ===================================================================
   float matrix2d_max (float **x, int m, int n) {
   float y=x[0][0];
   for (int j=0; j<m; j++) {
   for (int i=0; i<n; i++)
       if (y<x[j][i]) y=x[j][i];
   }
   return y;
}
// ===================================================================
   float matrix2d_amax (float **x, int m, int n) {
   float y=x[0][0];
   for (int j=0; j<m; j++) {
   for (int i=0; i<n; i++)
       if (y<fabs(x[j][i])) y=fabs(x[j][i]);
   }
   return y;
}
// ===================================================================
   float **matrix2d_product (float **a, float **b, int m, int n) {
   float **c=matrix2d_init_float (0.0,m,n);
   for (int i=0; i<m; i++) {
       for (int j=0; j<n; j++)
           c[i][j]=a[i][j]*b[i][j];
   }
   return c;
}
// ===================================================================
   void Matrix2dMul (float **y, float **a, float **b, int m, int n, int l) {
   for (int j=0; j<m; j++) {
   for (int i=0; i<l; i++) {
       y[j][i]=0.0;
       for (int k=0; k<n; k++)
           y[j][i] += a[j][k]*b[k][i];
   }
   }
}
// ===================================================================
   void M2dVecMul (float *y, float **a, float *b, int m, int n) {
   for (int j=0; j<m; j++) {
       y[j]=0.0;
       for (int k=0; k<n; k++)
           y[j] += a[j][k]*b[k];
   }
}
// ===================================================================
   void VecM2dMul (float *y, float *a, float **b, int m, int n) {
   for (int j=0; j<m; j++) {
       y[j]=0.0;
       for (int k=0; k<n; k++)
           y[j] += a[k]*b[k][j];
   }
}
// ===================================================================
   void MatAdd (float **y, float **a, float **b, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y[j][i]=a[j][i]+b[j][i];
}
// ===================================================================
   void MatSub (float **y, float **a, float **b, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y[j][i]=a[j][i]-b[j][i];
}
// ===================================================================
   void MatMul (float **y, float **a, float **b, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y[j][i]=a[j][i]*b[j][i];
}
// ===================================================================
   void MatMin (float **y, float **a, float **b, int m, int n) {
   for (int j=0; j<m; j++) {
       for (int i=0; i<n; i++) {
           if (a[j][i]<b[j][i])
              y[j][i]=a[j][i];
           else
              y[j][i]=b[j][i];
       }
   }
}
// ===================================================================
   void MatMax (float **y, float **a, float **b, int m, int n) {
   for (int j=0; j<m; j++) {
       for (int i=0; i<n; i++) {
           if (a[j][i]<b[j][i])
              y[j][i]=b[j][i];
           else
              y[j][i]=a[j][i];
       }
   }
}
// ===================================================================
   void MatSet (float **x, float a, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++) 
       x[j][i]=a;
}
// ===================================================================
   void MatScale (float **y, float **x, float a, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y[j][i]=x[j][i]*a;
}
// ===================================================================
   void MatNorm (float **y, int m, int n) {
   float p=matrix2d_amax(y,m,n);
   if (p>0.0) {
      for (int j=0; j<m; j++)
      for (int i=0; i<n; i++)
          y[j][i] /= p;
   }
}
// ===================================================================
   void MatInv (float **y, float **x, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++) {
       if (x[j][i]!=0.0) 
          y[j][i]=1.0/x[j][i];
       else
          y[j][i]=0.0;
   }
}
// ===================================================================
   void MatAbs (float **y, float **x, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y[j][i]=fabs(x[j][i]);
}
// ===================================================================
   float MatDot (float **a, float **b, int m, int n) {
   float y=0.0;
   for (int j=0; j<m; j++) {
   for (int i=0; i<n; i++)
       y += a[j][i]*b[j][i];
   }
   return y;
}
// ====== =============================================================
   float MatAbsMean (float **a, int m, int n) {
   float y=0.0;
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       y += fabs(a[j][i]);
   y = y/(m*n);
   return y;
}
// ===================================================================
   void MatAdd1 (float **z, float **x, float **y, float a, float b, int m, int n) {
   for (int j=0; j<m; j++)
   for (int i=0; i<n; i++)
       z[j][i]=a*x[j][i]+b*y[j][i];
}
// ===================================================================
   void MatCopy (float **y, float **x, int m, int n) {
   for (int j=0; j<m; j++)
       memcpy(y[j],x[j],n*sizeof(float));
}
// ===================================================================




// ===================================================================
  void Mat3dCopy (float ***y, float ***x, int l1, int l2, int m1, int m2, int n1, int n2) {
  
  int kk=0;
  for (int k=l1; k<=l2; k++) {
      int jj=0;
      for (int j=m1; j<=m2; j++) {
          int ii=0;
          for (int i=n1; i<=n2; i++)
              y[kk][jj][ii++]=x[k][j][i];
          jj += 1;
      }
      kk += 1;
  }
  return; 
}
// ===================================================================
