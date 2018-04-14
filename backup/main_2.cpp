#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "readpara.h"
#include "wyarray.h"
#include "kktecqfunc.h"
#include "solverlin.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    int maxLine = 6;
    string scoreFileName("data/aaa.txt");
    string targetFileName("data/aab.txt");

    //string fileName("data/usr2msg_score_0319.txt");
    //int nusr = 13282, nmsg = 200;
      
    //string fileName("data/usr2msg_score_0319_p10msg.txt");
    //int nusr = 2909, nmsg = 10;
      
    //int maxLine = 16070;
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
    KKTEcqFunc kktf(rdf.score, rdf.bottom, rdf.upper, rdf.nusr, rdf.nmsg);
    kktf.init();

    // define the solver
    SolverLin solin;

    // run the solver
    solin.InnerPoint (kktf);

    kktf.printResult("output/u2p_result.txt");

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    rdf.free();
    kktf.free();
}
