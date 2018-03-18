#include <iostream>
#include <stdlib.h>
#include "inout.h"
#include "assign.h"
#include "wyarray.h"
#include "inversion.h"
using namespace std;

int main (int argc, char **argv) {

    int nusr = 3;
    int nmsg = 2;
    int minMsgId = 50925;
    int ns = nusr * nmsg;
    int nm = ns + nusr;
    double *score = new double[ns];

    inout exaio;
    cout << " --- score --- " << endl;
    exaio.readtext_example (score, nusr, nmsg, minMsgId);
    exaio.printArray2d(score, nusr, nmsg);


    double lambda = 0.001;
    double gamma = 0.1;
    double *target = new double[nmsg];
    double *actual = new double[nmsg];

    // target flow for each msg
    for (int j = 0; j < nmsg; j ++) {
	target[j] = 1.0;
    }

    cout << " --- parameters --- " << endl;
    cout << " lambda = " << lambda << endl;
    cout << "  gamma = " <<  gamma << endl;

    double *dx = new double[nm];
    double *x0 = new double[nm];
    double *x1 = new double[nm];
    double *b0 = new double[nm];
    double *b1 = new double[nm];


    // initial value (assignment)
    for (int i = 0; i < nusr; i ++) {
        for (int j = 0; j < nmsg; j ++) {
	    x1[i * nmsg + j] = 1.0 / nmsg;
	}
    }
    for (int i = 0; i < nusr; i ++) {
	x1[ns + i] = 0.1;
    }

    // CG parameters
    int iprint = 0;
    int ncg = 100;
    double eps = 1e-06;
    double wns = 0.01;

    int iter_out = 0, outer = 50;
    int iter_inn = 0, inner = 10;
    double step = 0.001;        // search steplength
    double alpha = 0.01;
    double beta = 0.5;          // steplength shrinkage
    double eps0 = 1e-06;
    double r_perc = 1.0;
    double r_x0 = 1.0, r_x1, dr;
    double mis0, mis1;
    double yield;


    while (r_perc > eps0 && iter_out < outer) {

	iter_out += 1;
	cout << " ### " << endl;
	cout << " ### ------ outer=" << iter_out << "/" << outer << ":" << endl;

	memcpy(x0, x1, (nm) * sizeof(double));
	if (iprint == 1) {
	   cout << " --- initial x --- " << endl;
           exaio.printArray(x0, nm);
	}

	// calculate residual and norm
        assign ass(lambda, gamma, score, target, x0, nusr, nmsg);
        ass.residual(b0, x0);
        mis0 = ass.misfit(x0);
	r_x0 = norm2(b0, nm);

	if (iprint == 1) {
	   cout << " --- initial residual --- " << endl;
           exaio.printArray(b0, nm);
	}

        memset(dx, 0, (nm) * sizeof(double));

        solveLE sol(eps, wns, ncg, nm);
        sol.cgswy (dx, dx, b0, ass);
        //sol.cgsly (dx, dx, b0, ass);

	if (iprint == 1) {
	   cout << " --- initial dx --- " << endl;
           exaio.printArray(dx, nm);
	}

	iter_inn = 0;
	dr = 1.0;

	while (dr > 0 && iter_inn < inner) {
          iter_inn += 1;

	  /* x1 = x0 + step * dx; */
          VectorAdd (x1, x0, dx, 1.0, step, nm);
	  if (iprint == 1) {
	     cout << " --- update x ---- " << endl;
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
	        cout << " --- update residual ---- " << endl;
                exaio.printArray(b1, nm);
	     }
	  }

	  dr = r_x1 - (1.0 - alpha * step) * r_x0;
	  cout << " ------ inner=" << iter_inn << "/" << inner <<", step=" << step << ", mis=" << mis1 <<", r_x0=" << r_x0 << ", r_x1=" << r_x1 << ", dr=" << dr << endl;
          step = beta * step;

        }

	if (iter_out % 10 == 1) {
	   cout << " --- update x ---- " << endl;
           exaio.printArray(x1, nm);
	}

	// step = step / (beta * beta);
	
        ass.flow(actual, x1);
	yield=ass.yield(x1);
        r_perc = abs(r_x1 / r_x0 - 1.0);
	cout << " r perc = " << r_perc << endl;

    }

    cout << endl;
    cout << " --- final x ---- " << endl;
    exaio.printArray(x1, nm);
    cout << " --- flow x ---- " << endl;
    exaio.print2Array(actual, target, nmsg);
    cout << " yield=" << yield << endl;
    cout << " misfit=" << mis1 << endl;

    delete [] score;
    delete [] target;
    delete [] actual;
}
