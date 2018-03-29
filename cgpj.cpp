#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cgpj.h"
#include "wyarray.h"
#include "inout.h"
using namespace std;


/* ===================================================================
 * Projected CG method
 *
 * The original problem is
 *            min q(x) = (1/2) * x'Gx + x'c
 * s.t.             Ax = b
 *
 * KKT condition:
 * | G   A'| |-x|   | c|
 * |       |*|  | = |  |
 * | A   0 | | v|   |-b|
 * 
 * The matrix is indefinite, and CG is not appropriate for solving it.
 * the initial point x satisfies Ax=b
 * and r=Gx+c
 * ================================================================ */


//SolverPCG::SolverPCG (double eps, double wns, int niter, KKTEcqFunc *ass) {
//SolverPCG::SolverPCG (double eps, double wns, int niter, KKTMatrix *ass) {
SolverPCG::SolverPCG (double eps, double wns, int niter, KKTlp *ass) {
   this->eps = eps;
   this->wns = wns;
   this->niter = niter;
   this->n = ass->n;
   this->m = ass->m;
   this->ass = ass;
   this->iprint = 0;
}



/* solve v from AA'v = Ar 
 * g = r - A'v
 */
void SolverPCG::ProjOpr (double *g, double *r) {
    double *v = new double[m];
    double *p = new double[m];
    double *q = new double[n];
    memset(v, 0, m * sizeof(double));
    memset(p, 0, m * sizeof(double));
    memset(q, 0, n * sizeof(double));
    ass->Forward(p, r);
    if (iprint == 1) ass->PrintVector(p, m, "PCG ProjOpr p");

    this->SolverCg(v, p, 0);
    if (iprint == 1) ass->PrintVector(v, m, "PCG ProjOpr v");

    ass->Adjoint(q, v);
    VectorSub(g, r, q, n);
    delete [] p;
    delete [] q;
    delete [] v;
}

    

/* yc = A'v - Gx 
 * yb = Ax        */
void SolverPCG::KKTfunc(double *yc, double *yb, double *x, double *v) {
    double *p = new double[n];
    memset(p, 0, n * sizeof(double));
    memset(yb, 0, m * sizeof(double));
    ass->Forward(yb, x);
    ass->Adjoint(p, v);
    ass->G_opr(yc, x);
    VectorSub(yc, p, yc, n);
    delete [] p;
}

  
/* A'Ax = A'b */
void SolverPCG::GetInitVo (double *x, double *b) {
    double *p = new double[n];
    memset(p, 0, n * sizeof(double));
    memset(x, 0, n * sizeof(double));
    //this->CheckCgsOp(1);
    ass->Adjoint(p, b);
    this->SolverCg(x, p, 1);
    delete [] p;
}

/* r = Gx + c */
void SolverPCG::GetResi(double *r, double *x, double *c) {
    ass->G_opr(r, x); 
    for (int i = 0; i < n; i ++)
	r[i] += c[i];
}

/* A'v = r = Gx + c */
void SolverPCG::GetVfromX (double *v, double *x, double *c) {
   double * r = new double[n];
   double * p = new double[m];
   memset(p, 0, m * sizeof(double));
   memset(v, 0, m * sizeof(double));
   this->GetResi(r, x, c);
   ass->Forward(p, r);
   this->SolverCg(v, p, 0);
   delete [] r;
   delete [] p;
}

void SolverPCG::ProjCG (double *x, double *v, double *c, double *b) {

   int nm = ass->n;
   double alpha, beta;
   double r1, r2, b2;
   double *y = new double[nm];
   double *r = new double[nm];
   double *g = new double[nm];
   double *p = new double[nm];

   this->GetInitVo(x, b);
   if (iprint == 1) ass->PrintVector(x, nm, "PCG init x");

   this->GetResi(r, x, c);
   if (iprint == 1) ass->PrintVector(r, nm, "PCG resi r");

   this->ProjOpr(g, r);
   if (iprint == 1) ass->PrintVector(g, nm, "PCG ProjCG g");

   for (int k=0; k<nm; k++) {
       p[k] = -g[k];
   }

   r2=dot_product(g, g, nm);
   b2=dot_product(c, c, nm);

   for (int k=0; k<this->niter; k++) {

      //if (k%5==0) 
      //cout << " r2=" << r2 << " b2=" << b2 << endl;
      //printf(" *** iter=%d, eps=%f\n", k+1, sqrt(r2/b2));
      if (r2<b2*this->eps*this->eps) break;

      ass->G_opr(y, p);
      //ass->PrintVector(y, nm, "PCG ProjCG y");

      r1 = r2;
      alpha = r1 / dot_product(p,y,nm);
      if (iprint == 1) cout << "alpha=" << alpha << endl;

      for (int i=0; i<nm; i++) {
          x[i] = x[i] + alpha * p[i];
          r[i] = r[i] + alpha * y[i];
      }
      if (iprint == 1) ass->PrintVector(r, nm, "PCG ProjCG r");

      this->ProjOpr(g, r);
      if (iprint == 1) ass->PrintVector(g, nm, "PCG ProjCG g");

      r2 = dot_product(r, g, nm);
      if (iprint == 1) cout << "PCG ProjCG r1=" << r1 << " r2=" << r2 << endl;

      beta = r2 / r1;

      for (int i=0; i<nm; i++)
          p[i] = -g[i] + beta * p[i];
    
   }

   this->GetVfromX (v, x, c);

   delete [] y;
   delete [] r;
   delete [] p;
   delete [] g;
}
// ===================================================================



/* solve AA'x = b */
/* solve A'Ax = b */
void SolverPCG::SolverCg(double *x, double *b, int flag) {
   double eps = 1e-06;
   int niter = 100;
   int nt = 0;
   if (flag == 1) 
      nt = ass->n;
   else if (flag == 0)
      nt = ass->m;


   double alpha, beta;
   double r1, r2, b2;
   double *y = new double[nt];
   double *r = new double[nt];
   double *p = new double[nt];

   ass->Opr(y, x, flag);
   if (iprint == 1) {
      cout << "flag=" << flag << " n=" << nt << endl;
      ass->PrintVector(x, nt, "PCG CG x");
      ass->PrintVector(y, nt, "PCG CG y");
      ass->PrintVector(b, nt, "PCG CG b");
   }
   for (int k=0; k<nt; k++) {
       r[k]=b[k]-y[k];
       p[k]=r[k];
   }

   r2=dot_product(r,r,nt);
   b2=dot_product(b,b,nt);
 
   // printf("   *** CG starts \n");
   for (int k=0; k<niter; k++) {

      //if (k%5==0) 
      //printf(" *** iter=%d, eps=%f\n", k+1, sqrt(r2/b2));
      if (r2<b2*eps*eps) break;

      ass->Opr(y, p, flag);
      r1=r2;
      alpha=r1/dot_product(p,y,nt);
      for (int i=0; i<nt; i++) {
          x[i]=x[i]+alpha*p[i];
          r[i]=r[i]-alpha*y[i];
      }
      r2=dot_product(r,r,nt);
      beta=r2/r1;
      for (int i=0; i<nt; i++)
          p[i]=r[i]+beta*p[i];
   }

   delete [] y;
   delete [] r;
   delete [] p;
} 



void SolverPCG::CheckCgsOp (int flag) {

   int nt = 0;
   if (flag == 1) 
      nt = ass->n;
   else if (flag == 0)
      nt = ass->m;

   cout << " nt=" << nt << " " << ass->m << endl;

   double *x1 = (double *)malloc(nt*sizeof(double));
   double *y1 = (double *)malloc(nt*sizeof(double));
   double *x2 = (double *)malloc(nt*sizeof(double));
   double *y2 = (double *)malloc(nt*sizeof(double));

   srand((unsigned)time(NULL));
   for (int k=0; k<nt; k++) {
       x1[k] = rand()/(RAND_MAX+1.0);
       x2[k] = rand()/(RAND_MAX+1.0);
   }

   ass->Opr (y1,x1,flag);     /* y1 = A x1 */
   ass->Opr (y2,x2,flag);     /* y2 = A x2 */

   double z1=dot_product(x1,y2,nt);
   double z2=dot_product(x2,y1,nt);

   double c1=dot_product(x1,y1,nt);
   double c2=dot_product(x2,y2,nt);

   printf(" n=%d\n",nt);
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
