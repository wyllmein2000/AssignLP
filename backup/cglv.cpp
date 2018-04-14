#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cglv.h"
#include "wyarray.h"
#include "inout.h"
using namespace std;


SolverCG::SolverCG (double eps, double wns, int niter, int nm) {
   this->eps = eps;
   this->wns = wns;
   this->niter = niter;
   this->nm = nm;
}


void SolverCG::init () {
   this->niter = 100;
   this->eps = 1e-06;
   this->wns = 0.01;
}


void SolverCG::cgstabilize (double *y, double *x) {
   double ym=vector_amax(y,nm);
   double xm=vector_amax(x,nm);
   if (xm>0) {
      for (int k=0; k<nm; k++)
          y[k] += wns*ym/xm*x[k];
   }
}


// ===================================================================
// G*x=t
// G'Gx=G't
// Ax=b
// s0 as an initial guess
// t as known data
// 
// Input:  x0, b
// Output: x=s=1/v
// ===================================================================
void SolverCG::CgsWy (double *x, double *x0, double *b, PDMtrFunc &ass) {

   double alpha, beta;
   double r1, r2, b2;
   double *y = new double[nm];
   double *r = new double[nm];
   double *p = new double[nm];

   // obtain y=A(v0)*v0
   // printf("   *** \n");
   ass.op(y, x0);
   for (int k=0; k<nm; k++) {
       r[k]=b[k]-y[k];
       p[k]=r[k];
       x[k]=x0[k];
   }

   r2=dot_product(r,r,nm);
   b2=dot_product(b,b,nm);
   // printf("   *** %f %f\n",r2,b2);

   // printf("   *** CG starts \n");
   for (int k=0; k<this->niter; k++) {

      //if (k%5==0) 
      //printf(" *** iter=%d, eps=%f\n", k+1, sqrt(r2/b2));
      if (r2<b2*this->eps*this->eps) break;

      // obtain y=G'(p)*G(p)*p
      ass.op(y, p);

      r1=r2;
      alpha=r1/dot_product(p,y,nm);
      for (int i=0; i<nm; i++) {
          x[i]=x[i]+alpha*p[i];
          r[i]=r[i]-alpha*y[i];
      }
      r2=dot_product(r,r,nm);
      beta=r2/r1;
      for (int i=0; i<nm; i++)
          p[i]=r[i]+beta*p[i];
    
   }

   delete [] y;
   delete [] r;
   delete [] p;
} 
// ===================================================================
void SolverCG::cgsly (double *x, double *x0, double *b, PDMtrFunc ass) {

   double alpha, beta;
   double r1, r2, b2;
   double xm, pm, ym;

   double *y = new double[nm];
   double *r = new double[nm];
   double *p = new double[nm];

   memcpy(x,x0,nm*sizeof(double));

   // obtain y=A(v0)*v0
   printf("   *** \n");

   //op (y,x);
   ass.op(y, x);
   cgstabilize (y,x);

   char istr[3];
   //writedata("x00.bin",x,nm); 
   //writedata("y00.bin",y,nm); 

   //
   for (int k=0; k<nm; k++) {
       r[k]=b[k]-y[k];
       p[k]=r[k];
       // printf("%d %f %f %f\n",k,y[k],b[k],r[k]);
   }

   b2=dot_product(b,b,nm);
   r2=dot_product(r,r,nm);
   printf("   *** yes %f %f\n",r2,b2);
   if (b2<=1e-06) return;


   printf("   *** CG starts \n");
   for (int k=1; k<=this->niter; k++) {

//    if (k%50==0)
      printf(" *** iter=%d, eps=%f\n", k, sqrt(r2/b2));
      if (r2<b2*this->eps*this->eps) break;

      // obtain y=G'(p)*G(p)*p
      ass.op(y, p);
      cgstabilize (y,p);


      r1=r2;
      alpha=r1/dot_product(p,y,nm);
      for (int i=0; i<nm; i++) {
          x[i]=x[i]+alpha*p[i];
          //r[i]=r[i]-alpha*y[i];
      }

/*
      snprintf(istr,3,"%2.2d",k);
      char fb1[30]="x";
      char *fb2=strcat(strcat(fb1,istr),".bin");
      writedata(fb2,x,nm); 
      char fb3[30]="y";
      char *fb4=strcat(strcat(fb3,istr),".bin");
      writedata(fb4,y,nm); 
*/

//    if (k==1) printf("alpha=%f\n",alpha);
//    if (k==1) writedata("rra.bin",r,nm);
//    if (k==1) writedata("yya.bin",y,nm);

      //
      if (k%150==148) {
         ass.op(y, x);
         cgstabilize (y,x);
         VectorSub (r,b,y,nm);
      }
      else {
         for (int i=0; i<nm; i++) r[i]=r[i]-alpha*y[i];
         //r=vector_sub(r,vector_scale(y,alpha,nm),nm);
      }

      //if (k==1) writedata("rrb.bin",r,nm);

      //
      r2=dot_product(r,r,nm);
      beta=r2/r1;
      for (int i=0; i<nm; i++)
          p[i]=r[i]+beta*p[i];
    

//    char fb3[30]="r";
//    char *fb4=strcat(strcat(fb3,istr),".bin");
//    writedata(fb4,r,nm); 
   }

   delete [] y;
   delete [] r;
   delete [] p;
} 
// ===================================================================
void SolverCG::cgsmy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*)) {

   double alpha, beta;
   double r1, r2, b2;
   double *y = vector_init_double (0.0,nm);
   double *r = vector_init_double (0.0,nm);
   double *r0= vector_init_double (0.0,nm);
   double *u = vector_init_double (0.0,nm);

   double *p = vector_init_double (0.0,nm);
   double *q = vector_init_double (0.0,nm);

   double ratio = 1.0;
   double rho1 = 1.0;
   double rho2;
   
   // obtain y=A(v0)*v0
   memcpy(x,x0,nm*sizeof(double));
   op (y,x);
   cgstabilize (y,x);

   //
   for (int k=0; k<nm; k++) r0[k]=b[k]-y[k];
   memcpy(r,r0,nm*sizeof(double));

   b2=dot_product(b,b,nm);
   r2=dot_product(r,r,nm);
   printf("   *** %f %f\n",r2,b2);

// writedata("by.bin",b2,nm);
// writedata("ry.bin",r2,nm);

   //
   printf("   *** CG starts \n");
   for (int k=0; k<niter; k++) {

//    if (k%5==0)
         printf(" *** iter=%d, eps=%f, ratio=%f\n", k+1, sqrt(r2/b2), ratio);

      if (r2<b2*eps*eps || ratio<1e-08) break;

      //
      rho2=dot_product(r0,r,nm);
      beta=rho2/rho1;
      rho1=rho2;

//    printf("%d %f %f %f\n",k,beta,rho1,rho2);

      for (int i=0; i<nm; i++) {
          u[i]=r[i]+beta*q[i];
          p[i]=u[i]+beta*(q[i]+beta*p[i]);
      }


      //print1d(u,nm);
      //print1d(p,nm);

      // obtain y=G'(p)*G(p)*p
      op (y,p);
      cgstabilize (y,p);

      //
      alpha=rho2/dot_product(r0,y,nm);
      for (int i=0; i<nm; i++) {
          q[i]=u[i]-alpha*y[i];
          x[i]=x[i]+alpha*(u[i]+q[i]);
      }

      op (y,x);
      cgstabilize (y,x);
      VectorSub (r,b,y,nm);    // r=b-y;


      //
      r1=r2;
      r2=dot_product(r,r,nm);
      ratio=fabs(1.0-r2/r1);
   }

   free(r0);
   free(r);
   free(y);
   free(u);
   free(p);
   free(q);
} 
// ===================================================================
void SolverCG::bcgswy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*)) {

   double r2, b2;
   double rho,w,beta;

   memcpy(x,x0,nm*sizeof(double));

   printf("   *** \n");
   double *yy = vector_init_double (0.0,nm);
   double *y = vector_init_double (0.0,nm);
   double *p = vector_init_double (0.0,nm);
   double *v = vector_init_double (0.0,nm);
   double *s = vector_init_double (0.0,nm);

   op (y,x);
   // cgstabilize (y,x,wns,nm);
   double *r0 = (double *)malloc(nm*sizeof(double));
   double *r  = (double *)malloc(nm*sizeof(double));
   VectorSub (r0,b,y,nm);
   VectorSub (r, b,y,nm);
   //memcpy(r,r0,nm*sizeof(double));

   double rho0=1.0;
   double alpha=1.0;
   double w0=1.0;
   
   b2=dot_product(b,b,nm);
   r2=dot_product(r,r,nm);
   printf("   *** %f %f\n",r2,b2);

   //
   printf("   *** CG starts \n");
   for (int k=0; k<niter; k++) {

      if (k%50==0)
         printf(" *** iter=%d, eps=%f\n", k+1, sqrt(r2/b2));
      if (r2<b2*eps*eps) {
         printf(" *** %d iterations ...\n",k);
         break;
      }

      //
      rho=dot_product(r0,r,nm);
      beta=(rho/rho0)*alpha/w0;

      for (int i=0; i<nm; i++)
          p[i]=r[i]+beta*(p[i]-w0*v[i]);

      // obtain y=G'(p)*G(p)*p
      op (v,p);
      // cgstabilize (y,p,wns,nm);

      alpha=rho/dot_product(r0,v,nm);
      for (int i=0; i<nm; i++)
          s[i]=r[i]-alpha*y[i];

      op (y,s);
      w=dot_product(y,s,nm)/dot_product(y,y,nm);

      for (int i=0; i<nm; i++)
          x[i]=x[i]+alpha*p[i]+w*s[i];

      
      //
      op (yy,x);
      for (int i=0; i<nm; i++) 
          r[i]=b[i]-yy[i];
      r2=dot_product(r,r,nm);
      //

      for (int i=0; i<nm; i++) 
          r[i]=s[i]-w*y[i];
      r2=dot_product(r,r,nm);

      w0=w;
      rho0=rho;

   }

   free(r0);
   free(r);
   free(y);
   free(p);
   free(v);
   free(s);
   free(yy);
} 


// ===================================================================
void SolverCG::CheckCgsOp (PDMtrFunc &ass) {

   int n = this->nm;
   double *x1 = (double *)malloc(n*sizeof(double));
   double *y1 = (double *)malloc(n*sizeof(double));
   double *x2 = (double *)malloc(n*sizeof(double));
   double *y2 = (double *)malloc(n*sizeof(double));

   srand((unsigned)time(NULL));
   for (int k=0; k<n; k++) {
       x1[k] = rand()/(RAND_MAX+1.0);
       x2[k] = rand()/(RAND_MAX+1.0);
   }

// printf("hello\n");
   ass.op (y1,x1);     /* y1 = A x1 */
   ass.op (y2,x2);     /* y2 = A x2 */

// printf("hellu\n");
   double z1=dot_product(x1,y2,n);
   double z2=dot_product(x2,y1,n);

   double c1=dot_product(x1,y1,n);
   double c2=dot_product(x2,y2,n);

   printf(" n=%d\n",n);
   printf(" *** check symmetry *** \n");
   printf("x1' * A * x2 = %f\n",z1); 
   printf("x2' * A * x1 = %f\n",z2); 
   printf("\n");
  
   printf(" *** check positive *** \n");
   printf("x1' * A * x1 = %f\n",c1); 
   printf("x2' * A * x2 = %f\n",c2); 
   printf("\n");

   free(x1);
   free(y1);
   free(x2);
   free(y2);
}
// ===================================================================
/*
   void checksym (int n, double (*op)(double*, double*)) {

   double *x1 = (double *)malloc(n*sizeof(double));
   double *y1 = (double *)malloc(n*sizeof(double));
   double *x2 = (double *)malloc(n*sizeof(double));
   double *y2 = (double *)malloc(n*sizeof(double));

   srand((unsigned)time(NULL));
   for (int k=0; k<n; k++) {
       x1[k] = rand()/(RAND_MAX+1.0);
       x2[k] = rand()/(RAND_MAX+1.0);
   }

   op (y1,x1);
   op (y2,x2);


   double z1=dot_product(x1,y2,n);
   double z2=dot_product(x2,y1,n);

   printf("y2'*x1= %f\n",z1); 
   printf("y1'*x2= %f\n",z2); 
  
   free(x1);
   free(y1);
   free(x2);
   free(y2);
}*/
// ===================================================================


















// ===================================================================
   int m=4;
   int n=5;
   double **g=matrix2d_init_double(0.0,m,n);
// ===================================================================
// y=g'x
// ===================================================================
void SolverCG::matadj (double *y, double *x) {
   for (int i=0; i<n; i++) {
       y[i]=0.0;
       for (int j=0; j<m; j++)
           y[i] += g[j][i]*x[j];
   }
}
// ===================================================================
void SolverCG::matop (double *y, double *x) {
   double *z=vector_init_double(0.0,m);
   M2dVecMul (z,g,x,m,n);
   matadj (y,z);
   free(z);
}
// ===================================================================
void SolverCG::checkcgs () {

   int niter=10;
   double eps=0.01;
   double wns=0.01;

   double *x=vector_init_double(0.0,n);
   double *t=vector_init_double(0.0,m);

   for (int i=0; i<n; i++) {
       for (int j=0; j<m; j++)
           g[j][i]=2*j+i-5.0;
       x[i]=i; 
   }

   M2dVecMul(t,g,x,m,n);

   double *x0=vector_init_double(0.0,n);
   double xa=vector_sum(x,n)/n;
   for (int i=0; i<n; i++) x0[i]=xa;

   double *b=vector_init_double(0.0,n);
   matadj (b,t);

   double *xx=vector_init_double(0.0,n);

// ---
   //cgswy (xx,x0,b,eps,niter,n,matop);
   printf("  cgs wuyan: \n");
   for (int i=0; i<n; i++)
       printf("i=%d,%f,%f\n",i,x[i],xx[i]);

// ---
   printf("  cgs luoyi: \n");
   memset(xx,0,n*sizeof(double));
   //cgsly (xx,x0,b,wns,eps,niter,n,matop);
   for (int i=0; i<n; i++)
       printf("i=%d,%f,%f\n",i,x[i],xx[i]);

// ---
   printf("  cgs mayue: \n");
   memset(xx,0,n*sizeof(double));
   //cgsmy (xx,x0,b,wns,eps,niter,n,matop);
   for (int i=0; i<n; i++)
       printf("i=%d,%f,%f\n",i,x[i],xx[i]);

  
// ---
   printf("  bcgs wuyan: \n");
   memset(xx,0,n*sizeof(double));
   //bcgswy (xx,x0,b,wns,eps,niter,n,matop);
   for (int i=0; i<n; i++)
       printf("i=%d,%f,%f\n",i,x[i],xx[i]);


}
// ===================================================================
// ===================================================================
  
