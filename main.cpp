#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "assign.h"
#include "wyarray.h"
#include "cglv.h"
using namespace std;

int main (int argc, char **argv) {
   
    double start = time(NULL);
    string fileName("data/aaa.txt");
    int nusr = 6, nmsg = 2;

    //string fileName("data/usr2msg_score_0319.txt");
    //int nusr = 13282, nmsg = 200;
    //
    //string fileName("data/usr2msg_score_0319_p10msg.txt");
    //int nusr = 2909, nmsg = 10;
    //
    //string fileName("data/u2p_score_0319.txt");
    //int nusr = 1000, nmsg = 17;

    int ns = nusr * nmsg;
    int nm = ns + nusr;
    double *score = new double[ns];
    memset(score, 0, ns * sizeof(double));

    inout exaio(fileName, nusr, nmsg);
    cout << " --- score --- " << endl;
    //exaio.readtext_example (score, nusr, nmsg, minMsgId);
    //exaio.printArray2d(score, nusr, nmsg);

    exaio.readtext_usr_msg_score (score);
    //exaio.printArray2d(score, nusr, nmsg);
    //exit(0);

    double lambda = 1e-03;
    double gamma = 1e-03;
    double *target = new double[nmsg];
    double *actual = new double[nmsg];
    double *actual_r = new double[nmsg];

    // target flow for each msg
    // srand(seed);
    for (int j = 0; j < nmsg; j ++) {
	//target[j] = 60.0 + (rand() % 101 - 50);
	target[j] = 3.0;
    }
    //target[0] = 4; target[1] = 4.0;

    cout << " --- parameters --- " << endl;
    cout << " lambda = " << lambda << endl;
    cout << "  gamma = " <<  gamma << endl;

    double *dx = new double[nm];
    double *x0 = new double[nm];
    double *x1 = new double[nm];
    double *xr = new double[nm];
    double *b0 = new double[nm];
    double *b1 = new double[nm];


    // initial value (assignment)
    for (int i = 0; i < nusr; i ++) {
        for (int j = 0; j < nmsg; j ++) {
	    x1[i * nmsg + j] = 1.0 / nmsg;
	    //x1[i * nmsg + j] = 0.1 * (1 + rand() % 9);
	}
    }
    for (int i = 0; i < nusr; i ++) {
	x1[ns + i] = 0.001;
    }

    // CG parameters
    int iprint = 0;
    int ncg = 100;
    double eps = 1e-06;
    double wns = 0.01;

    int iter_out = 0, outer = 50;
    int iter_inn = 0, inner = 10;
    double step = 1e-03;        // search steplength
    double alpha = 0.01;
    double beta = 0.5;          // steplength shrinkage
    double eps0 = 1e-06;
    double r_perc = 1.0;
    double r_x0 = 1.0, r_x1, dr;

    double mis0, mis1;
    double entropy, yield;
    double entropy_r, yield_r;
    double r1, r2;

    assign ass(lambda, gamma, score, target, x0, nusr, nmsg);
    solveLE sol(eps, wns, ncg, nm);

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
        ass.residual(b0, x0);
        mis0 = ass.misfit(x0);
	r_x0 = norm2(b0, nm);

	if (iprint == 1) {
	   cout << " --- initial residual --- " << endl;
           exaio.printArray(b0, nm);
	}

        memset(dx, 0, (nm) * sizeof(double));

        sol.cgswy (dx, dx, b0, ass);
        //sol.cgsly (dx, dx, b0, ass);

	if (iprint == 1) {
	   cout << " --- initial dx --- " << endl;
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

	if (iter_out % 10 == 1 && iprint == 1) {
	   cout << " --- update x ---- " << endl;
           exaio.printArray(x1, nm);
	}

	// step = step / (beta * beta);
	
        r_perc = abs(r_x1 / r_x0 - 1.0);
	cout << " --- outer = " << iter_out << "/" << outer << " mis=" << mis1 << " r perc = " << r_perc << endl;

	ass.maxres(&r1, &r2, b1);
        ass.flow(actual, x1);
	entropy = ass.entropy(actual, nmsg);
	yield = ass.yield(x1);

	ass.round(xr, x1);
        ass.flow(actual_r, xr);
	entropy_r = ass.entropy(actual_r, nmsg);
	yield_r = ass.yield(xr);

    }

    cout << endl;
    cout << " --- final x ---- " << endl;
    //exaio.printArray(x1, nm);
    cout << endl << " --- original flow x ---- " << endl;
    exaio.print2Array(actual, target, nmsg);
    cout << " res_1 = " << r1 << " res_2 = " << r2 << endl;
    cout << " yield=" << yield <<" entropy=" << entropy << " misfit=" << mis1 << endl;

    cout << endl << " --- round() flow x ---- " << endl;
    exaio.print2Array(actual_r, target, nmsg);
    cout << " yield=" << yield_r <<" entropy=" << entropy_r << endl;

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    exaio.printResult(yield, yield_r, entropy, entropy_r, actual, actual_r, target, nmsg);

    delete [] score;
    delete [] target;
    delete [] actual;
    delete [] actual_r;
}
