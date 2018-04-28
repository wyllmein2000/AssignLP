#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "readpara.h"
#include "wyarray.h"
#include "solver_assign_dual.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    //int maxLine = 6;
    //string scoreFileName("data/aaa.txt");
    //string targetFileName("data/aab.txt");

    //int maxLine = 37962;     // 17 msgs, 10000 users
    //string scoreFileName("data/u2p_score_0412.txt");
    //string targetFileName("data/u2p_target_0412.txt");
      
    //string fileName("data/usr2msg_score_0319_p10msg.txt");
      
    /*
    int maxLine = 16070;   // 17 msgs, 1000 users
    //int maxLine = 1991;   // 2 msgs, 1000 users
    string scoreFileName("data/u2p_score_0319.txt");
    string targetFileName("data/u2p_target_0319.txt");
    string logFilename("output/u2p_dual0319_log.txt");
    string outFilename("output/u2p_dual0319_result.txt");
    */

    //
    int maxLine = 5089149;  // 391 msgs, 370498 users
    string scoreFileName("data/u2p_score_0424.txt");
    string targetFileName("data/u2p_target_0424.txt");
    string logFilename("output/u2p_0424_dual_0_log.txt");
    string outFilename("output/u2p_0424_dual_0_result.txt");
    //

    /* read parameters from file */
    ReaderFileAss *rdf = new ReaderFileAss(scoreFileName, targetFileName, maxLine); 
    rdf->getDataDimension ();
    rdf->readScoreFromFile ();
    //rdf->readTargetFromFile ();
    rdf->setTarget (0);
    //rdf->printScore();
    rdf->printTarget();
    //exit(0);

    int nusr = rdf->nusr;
    int nmsg = rdf->nmsg;
    int nx = nusr * nmsg;

    // define varibles
    double *x = new double[nx];

    // define the solver
    SolverAssDual *sol = new SolverAssDual(rdf->score, rdf->bottom, nusr, nmsg);
    sol->SetFilename(outFilename, logFilename);

    // run the solver
    sol->GradientDescent (x);

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    if (durationTime < 60)
       cout << endl << "Cost：" << durationTime << " sec" << endl;
    else if (durationTime < 3600)
       cout << endl << "Cost：" << durationTime/60.0 << " min" << endl;
    else
       cout << endl << "Cost：" << durationTime/3600.0 << " hour" << endl;

    ofstream fp(outFilename, ios::app);
    if (durationTime < 60)
       fp << " Cost : " << durationTime << " sec" << endl;
    else if (durationTime < 3600)
       fp << "Cost：" << durationTime/60.0 << " min" << endl;
    else
       fp << "Cost：" << durationTime/3600.0 << " hour" << endl;
    fp.close();

    delete [] x;
    rdf->free();
}
