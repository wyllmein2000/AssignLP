#ifndef CGLV_H
#define CGLV_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "pdmtrfunc.h"


class SolverCG {

   private:
       int nm;
       int niter;
       double eps;
       double wns;

   public:
       SolverCG (double eps, double wns, int niter, int nm);

       void init ();

       void cgstabilize (double *y, double *x);
       void CgsWy (double *x, double *x0, double *b, PDMtrFunc &ass);
       void cgsly (double *x, double *x0, double *b, PDMtrFunc ass);

       void cgsmy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void bcgswy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void CheckCgsOp (PDMtrFunc &ass);

       void matadj (double *y, double *x);
       void matop (double *y, double *x);
       void checkcgs ();
}
;

#endif
