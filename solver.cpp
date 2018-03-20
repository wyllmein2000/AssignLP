#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "assign.h"
#include "wyarray.h"
#include "solver.h"
using namespace std;


/* -------------------------------------------------
 * algorithm
 *
 * while (r2/r1 - 1 > eps)
 *    compute dx
 *    while (dr > 0)
 *       step = step * beta
 *       x = x0 + step * dx
 *       dr = r2 - (1 - alpha * step) * r1
   ---------------------------------------------- */
SolverNle::SolverNle (int nm) {
    this->nm = nm;
    this->inner = 10;
    this->outer = 50;

    // x1 = x0 + step * dx
    this->step = 0.001;
    this->alpha = 0.1;
    this->beta = 0.5;
    this->eps = 1e-06;

    this->iprint = 0;

}


// equality constrained minimization
void SolverNle::NewtonEcm (assign &ass) {

    double *dx = new double[nm];
    double *x0 = new double[nm];
    double *x1 = new double[nm];

    double *xr = new double[nm];
    double *b0 = new double[nm];
    double *b1 = new double[nm];


    // CG parameters
    int cg_niter = 50;
    double cg_eps = 1e-06;
    double cg_wns = 0.01;

    // define a CG operator used in computing dx
    SolverCG solcg (cg_eps, cg_wns, cg_niter, nm);

    // check the operator
    // solcg.CheckCgsOp(ass); exit(0);


    // define an inout class
    inout exaio;

    // copy initial value from ass to x1
    ass.GetX (x1);

    int iter_out = 0;
    int iter_inn = 0;
    double r_perc = 1.0;
    double r_x0 = 1.0, r_x1, dr;

    double mis0, mis1;
    double r1, r2;

    while (r_perc > this->eps && iter_out < this->outer) {

	iter_out += 1;
	cout << " ### " << endl;
	cout << " ### ------ outer=" << iter_out << "/" << outer << ":" << endl;
        memcpy(x0, x1, nm * sizeof(double));

	if (iprint == 1) {
	   cout << endl << " --- initial x --- " << endl;
           exaio.printArray(x0, nm);
	}

	// calculate residual and norm
        ass.residual(b0, x0);
        mis0 = ass.misfit(x0);
	r_x0 = norm2(b0, nm);



	if (iprint == 1) {
	   cout << endl << " --- initial residual --- " << endl;
           exaio.printArray(b0, nm); exit(0);
	}



        memset(dx, 0, nm * sizeof(double));

        solcg.CgsWy (dx, dx, b0, ass);
        //solcg.cgsly (dx, dx, b0, ass);

	if (iprint == 1) {
	   cout << endl << " --- solve dx --- " << endl;
           exaio.printArray(dx, nm);
	}

	iter_inn = 0;
	dr = 1.0;
	step /= beta;

	while (dr > 0 && iter_inn < inner) {
          iter_inn += 1;

	  /* x1 = x0 + step * dx; */
          VectorAdd (x1, x0, dx, 1.0, step, nm);
	  if (iprint == 1) {
	     cout << endl << " --- inner update x ---- " << endl;
             exaio.printArray(x1, nm);
	  }

	  if (vector_min(x1, nm) < 1e-36) {
             mis1 = 0.0;
             r_x1 = 2.0 * r_x0;
	  }
          else {
             ass.residual(b1, x1);
             mis1 = ass.misfit(x1);
	     r_x1 = norm2(b1, nm);
	     if (iprint == 1) {
	        cout << endl << " --- update residual ---- " << endl;
                exaio.printArray(b1, nm);
	     }
	  }

	  dr = r_x1 - (1.0 - alpha * step) * r_x0;
	  cout << " ------ inner=" << iter_inn << "/" << inner <<", step=" << step << ", mis=" << mis1 <<", r_x0=" << r_x0 << ", r_x1=" << r_x1 << ", dr=" << dr << endl;
          step = beta * step;

        }

	if (iter_out % 10 == 1 && iprint == 1) {
	   cout << " --- outer update x ---- " << endl;
           exaio.printArray(x1, nm);
	}

	// update ass.x0
	ass.UpdateX (x1);
	ass.maxres(&r1, &r2, b1);
        r_perc = abs(r_x1 / r_x0 - 1.0);
	cout << " --- outer = " << iter_out << "/" << outer << " mis=" << mis1 << " r_perc = " << r_perc << " r1=" << r1 << " r2=" << r2 << endl;


    }

    cout << endl;
    cout << " --- final x ---- " << endl;
    //exaio.printArray(x1, nm);
    
    delete [] dx;
    delete [] x0;
    delete [] x1;
    delete [] xr;
    delete [] b0;
    delete [] b1;
}
