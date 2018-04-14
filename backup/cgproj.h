#ifndef INVERSION_H
#define INVERSION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "kktecqfunc.h"


class SolverPCG {

   private:
       int nm;
       int niter;
       double eps;
       double wns;

   public:
       SolverPCG (double eps, double wns, int niter, int nm);

       void Init ();
       void ProjCG (double *x, double *x0, double *b, KKTEcqFunc &ass);
       //void CheckCgsOp (assign &ass);

};

#endif
