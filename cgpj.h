#ifndef CGPJ_H
#define CGPJ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include "kktecqfunc.h"
#include "kktlp.h"
#include "kktmatrix.h"


class SolverPCG {

   private:
       int m, n;
       int niter;
       int iprint;
       double eps;
       double wns;
       //KKTEcqFunc *ass;
       //KKTMatrix *ass;
       KKTlp *ass;

   public:
       //SolverPCG (double eps, double wns, int niter, KKTEcqFunc *ass);
       //SolverPCG (double eps, double wns, int niter, KKTMatrix *ass);
       SolverPCG (double eps, double wns, int niter, KKTlp *ass);

       void Init ();

       void ProjOpr (double *, double *);
       void KKTfunc(double *, double *, double *x, double *v);
       void GetInitVo (double *, double *);
       void GetResi (double *, double *, double *);
       void GetVfromX (double *, double *, double *);

       void ProjCG (double *, double *, double *, double *);
       void SolverCg (double *, double *, int);
       void CheckCgsOp (int);
       //void CheckCgsOp (assign &ass);

};
#endif
