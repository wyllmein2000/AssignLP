#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "kktlp.h"
#include "cgpj.h"
#include "wyarray.h"
using namespace std;


class SolverLin {
    private:
	double mu, eta;
        int niter;
        int inner;
	int outer;
	int iprint;

    public:
	int m, n;
	KKTlp *ass;

        SolverLin (KKTlp *);
	void GetInitPoint(double *, double *, double *, double *, double *, SolverPCG *);
	void ModifyFuncPara(double *, double *, double *);
	void UpdateDs(double *, double *);
	void GetRhsPredictor(double *, double *, double *, double *);
	void GetRhsCorrector(double *, double *, double *, double *, double *);
	void GetRhsBC(double *, double *, double *, double *, double *);
	void EstimateDr(double *, double *, double *);
	void UpdateXVS(double *, double *, double *, double *, double *, double *);
	void InnerPoint (double *, double *, double *, double *, double *);

        void TestPCG(SolverPCG *sol);
};
