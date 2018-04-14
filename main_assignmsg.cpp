#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "readpara.h"
#include "wyarray.h"
#include "proassign.h"
#include "solvernle.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    //int maxLine = 6;
    //string scoreFileName("data/aaa.txt");
    //string targetFileName("data/aab.txt");

    int maxLine = 37962;
    string scoreFileName("data/u2p_score_0412.txt");
    string targetFileName("data/u2p_target_0412.txt");
      
    //string fileName("data/usr2msg_score_0319_p10msg.txt");
      
    //int maxLine = 16070;   // 17 msgs
    //int maxLine = 1991;   // 2 msgs
    //string scoreFileName("data/u2p_score_0319.txt");
    //string targetFileName("data/u2p_target_0319.txt");


    /* read parameters from file */
    ReaderFileAss *rdf = new ReaderFileAss(scoreFileName, targetFileName, maxLine); 
    rdf->getDataDimension ();
    rdf->readScoreFromFile ();
    rdf->readTargetFromFile ();
    //rdf->printScore();
    //rdf->printTarget();
    //exit(0);

    int nx = rdf->nusr * rdf->nmsg;

    // define the specific problem
    ProAssign *ass = new ProAssign(rdf->score, rdf->bottom, rdf->nusr, rdf->nmsg);

    // define varibles
    double *x = new double[nx];
    double *x0 = new double[nx];
    double *ce = new double[rdf->nusr];
    double *cn = new double[rdf->nmsg];

    ass->GetInitPoint(x0, ce, cn);
    memcpy(x, x0, nx * sizeof(double));

    // define the solver
    SolverNle *sol = new SolverNle(ass, "output/msg_output.txt");

    // run the solver
    sol->AugLag (x, ce, cn);

    ass->printResult("output/msg_result.txt", x0, x);

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    rdf->free();
}
