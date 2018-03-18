#ifndef INVERSION_H
#define INVERSION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "assign.h"


class solveLE {

   private:
       int nm;
       int niter;
       double eps;
       double wns;

   public:
       solveLE (double eps, double wns, int niter, int nm);

       void cgstabilize (double *y, double *x);
       void cgswy (double *x, double *x0, double *b, assign);
       void cgsly (double *x, double *x0, double *b, assign);

       void cgsmy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void bcgswy (double *x, double *x0, double *b, double wns, double eps, int niter, int nm, void (*op)(double*, double*));

       void checkcgsop (int n, void (*op)(double*, double*));

       void matadj (double *y, double *x);
       void matop (double *y, double *x);
       void checkcgs ();
}
;

#endif
