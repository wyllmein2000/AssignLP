#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "wyarray.h"
#include "readpara.h"
#include "pdmatrix.h"
#include "cglv.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    int maxLine = 6;
    string scoreFileName("data/aaa.txt");
    string targetFileName("data/aab.txt");

    //int maxLine = 1991;   // two msgs
    //string scoreFileName("data/u2p_score_0319.txt");
    //string targetFileName("data/u2p_target_0319.txt");


    /* read parameters from file */
    ReaderFileAss rdf(scoreFileName, targetFileName, maxLine); 
    rdf.getDataDimension ();
    rdf.readScoreFromFile ();
    rdf.readTargetFromFile ();
    rdf.printScore();
    rdf.printTarget();
    //exit(0);


    // define the specific problem
    PDMatrix *pdm = new PDMatrix(rdf.score, rdf.nusr, rdf.nmsg);

    // set a solution
    double *x0 = new double[rdf.nmsg];
    double *x = new double[rdf.nmsg];
    double *b = new double[rdf.nmsg];
	    
    // define the solver
    SolverCG solcg(0.001, 0.01, 50, pdm);
    solcg.CheckCgsOp ();

    // 
    x[0] = 0.123;
    x[1] = 0.498;
    pdm->Opr(b, x);
    pdm->PrintVector(b, rdf.nmsg, "bbb");

    // solve the problem
    memset(x0, 0, rdf.nmsg * sizeof(double));
    solcg.CgsWy(x0, x0, b);

    cout << "index   predicted   true " << endl;
    for (int i = 0; i < rdf.nmsg; i ++) {
        cout << " i = " << i << ", " << x0[i] << ", " << x[i] << endl; 
    }

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    pdm->Free();
}
