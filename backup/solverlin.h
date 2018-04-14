#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "kktecqfunc.h"
#include "wyarray.h"
#include "cglv.h"
#include "cgproj.h"
using namespace std;


class SolverLin {
    private:
	int nm;
        int inner;
	int outer;
	int iprint;
        double step;
        double alpha;
        double beta;
        double eps;

    public:
        SolverLin ();
	void InnerPoint (KKTEcqFunc &ass);
};
