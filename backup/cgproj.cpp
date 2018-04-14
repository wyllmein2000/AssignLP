#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cgproj.h"
#include "wyarray.h"
#include "inout.h"
using namespace std;


SolverPCG::SolverPCG (double eps, double wns, int niter, int nm) {
   this->eps = eps;
   this->wns = wns;
   this->niter = niter;
   this->nm = nm;
}


void SolverPCG::Init () {
   this->niter = 100;
   this->eps = 1e-06;
   this->wns = 0.01;
}



/* ===================================================================
 * Projected CG method
 *
 * The original problem is
 *            min q(x) = (1/2) * x'Gx + x'c
 * s.t.             Ax = b
 *
 * KKT condition:
 *
 * | G   A'| |-dx|   |Gx + c|   |r|
 * |       |*|   | = |      | = | |
 * | A   0 | | v |   |Ax - b|   |h|
 * 
 * x(k+1) = x(k) + dx
 * 
 * The matrix is indefinite, and CG is not appropriate for solving it.
 *
 * the initial point x0 satisfies Ax=b
// ================================================================ */
void SolverPCG::ProjCG (double *x, double *x0, double *b, KKTEcqFunc &ass) {

   double alpha, beta;
   double r1, r2, b2;
   double *y = new double[nm];
   double *r = new double[nm];
   double *g = new double[nm];
   double *p = new double[nm];

   for (int k=0; k<nm; k++) {
       x[k] = x0[k];
       r[k] = b[k];
   }
   ass.P_opr(g, r);
   for (int k=0; k<nm; k++) {
       p[k] = -g[k];
   }

   r2=dot_product(g, g, nm);
   b2=dot_product(b, b, nm);
   cout << " r2=" << r2 << " b2=" << b2 << endl;

   for (int k=0; k<this->niter; k++) {

      //if (k%5==0) 
      printf(" *** iter=%d, eps=%f\n", k+1, sqrt(r2/b2));
      if (r2<b2*this->eps*this->eps) break;

      ass.G_opr(y, p);

      r1 = r2;
      alpha = r1 / dot_product(p,y,nm);

      for (int i=0; i<nm; i++) {
          x[i] = x[i] + alpha * p[i];
          r[i] = r[i] + alpha * y[i];
      }

      ass.P_opr(g, r);

      r2 = dot_product(r, g, nm);

      beta = r2 / r1;

      for (int i=0; i<nm; i++)
          p[i] = -g[i] + beta * p[i];
    
   }

   delete [] y;
   delete [] r;
   delete [] p;
}
// ===================================================================
















  
