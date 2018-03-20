#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "assign.h"
#include "wyarray.h"
#include "cglv.h"
using namespace std;


class SolverNle {
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
        SolverNle (int nm);
        void NewtonEcm (assign &ass);
};
