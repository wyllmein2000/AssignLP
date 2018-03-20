#ifndef INVERSION_H
#define INVERSION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "assign.h"


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
       void CgsWy (double *x, double *x0, double *b, assign &ass);
       void cgsly (double *x, double *x0, double *b, assign ass);

       void cgsmy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void bcgswy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void CheckCgsOp (assign &ass);

       void matadj (double *y, double *x);
       void matop (double *y, double *x);
       void checkcgs ();
}
;

#endif
