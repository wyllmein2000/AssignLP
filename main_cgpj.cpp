#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "wyarray.h"
#include "readpara.h"
#include "kktmatrix.h"
#include "cgpj.h"
using namespace std;


int main (int argc, char **argv) {
   
    double start = time(NULL);
    srand((unsigned)time(NULL));

    int n = 3;
    int m = 2;
    int ng = n * n;
    int na = m * n;

    double *G = new double[ng];   // size(G) = n x n
    double *A = new double[na];   // size(A) = m x n

    for (int i = 0; i < ng; i ++)
	G[i] = rand()/(RAND_MAX + 1.0);
    for (int i = 0; i < na; i ++)
	A[i] = rand()/(RAND_MAX + 1.0);
        
    // define the specific problem
    KKTMatrix *kktm = new KKTMatrix(G, A, m, n);
    kktm->PrintMatrix(G, n, n, "Gopr");
    kktm->PrintMatrix(A, m, n, "Aopr");

    // set a solution
    double *x = new double[n];
    double *v = new double[m];
    for (int i = 0; i < n; i ++)
        x[i] = rand()/(RAND_MAX + 1.0);
    for (int i = 0; i < m; i ++)
        v[i] = rand()/(RAND_MAX + 1.0);

    // get RHS of equation
    double *yc = new double[n];
    double *yb = new double[m];
    kktm->KKTfunc(yc, yb, x, v);
    kktm->PrintVector(yc, n, "ccc");
    kktm->PrintVector(yb, m, "bbb");
	    

    // define the solver
    SolverPCG solpcg(0.001, 0.01, 50, kktm);
    //solpcg.CheckCgsOp (0);
    //solpcg.CheckCgsOp (1); exit(0);


    // solve the problem
    double *x0 = new double[n];
    double *v0 = new double[m];
    solpcg.ProjCG(x0, v0, yc, yb);

    // print result
    cout << "index   predicted   true " << endl;
    for (int i = 0; i < n; i ++) {
        cout << " i, x = " << i << ", " << x0[i] << ", " << x[i] << endl; 
    }
    for (int i = 0; i < m; i ++) {
        cout << " i, v = " << i << ", " << v0[i] << ", " << v[i] << endl; 
    }

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    kktm->Free();
}
