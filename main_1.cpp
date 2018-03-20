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

    // define the specific problem
    assign ass(fileName, nusr, nmsg);
    ass.init();

    // define the solver
    SolverNle solne (nm);

    // run the solver
    solne.NewtonEcm (ass);

    ass.printResult("output/a.txt");

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    ass.free();
}
