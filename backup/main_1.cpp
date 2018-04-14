#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "assign.h"
#include "wyarray.h"
#include "solver.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    //string fileName("data/aaa.txt");
    //int maxLine = 12;

    //string fileName("data/usr2msg_score_0319.txt");
    //int nusr = 13282, nmsg = 200;
      
    //string fileName("data/usr2msg_score_0319_p10msg.txt");
    //int nusr = 2909, nmsg = 10;
      
      
    string fileName("data/u2p_score_0319.txt");
    //int maxLine = 1991;   // two msgs
    int maxLine = 16070;

    // define the specific problem
    assign ass(fileName, maxLine);
    ass.init();

    // define the solver
    SolverNle solne;

    // run the solver
    solne.NewtonEcm (ass);

    ass.printResult("output/u2p_result.txt");

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    ass.free();
}
