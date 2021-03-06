#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "inout.h"
#include "wyarray.h"
#include "readpara.h"
#include "kktlp.h"
#include "pdmatrix.h"
#include "solverlin.h"
#include "cglv.h"
using namespace std;

void SetPara (double *A, double *b, double *c, int na) {
    // define the operator
    for (int i = 0; i < na; i ++)
	A[i] = rand()/(RAND_MAX + 1.0);

    A[0] = -0.973528;
    A[1] =  0.0919762;
    A[2] =  0.843787;
    A[3] = -0.52274;
    A[4] =  0.686307;
    A[5] = -0.762521;
    b[0] =  0.758965;
    b[1] =  0.648628;
    c[0] = -0.732288;
    c[1] =  1.06376;
    c[2] = -0.742943;
}

void SetSolution (double *x, double *v, double *s, double *A, double *b, double *c, int n, int m) {

    double *AT = new double[m * n];
    double *ct = new double[n];
    double *cc = new double[m];
    double *bb = new double[n];

    for (int i = 0; i < n; i ++) {
        x[i] = rand()/(RAND_MAX + 1.0);
        s[i] = rand()/(RAND_MAX + 1.0);
    }
    for (int i = 0; i < m; i ++)
        v[i] = rand()/(RAND_MAX + 1.0);

    PDMatrix *pdm = new PDMatrix(A, m, n);
    pdm->Adjoint(bb, b);
    SolverCG *sol = new SolverCG(0.001, 0.01, 50, pdm);
    sol->CheckCgsOp();
    sol->CgsWy(x, x, bb);

    for (int i = 0; i < m; i ++)
    for (int j = 0; j < n; j ++)
        AT[j * m + i] = A[i * n + j];
    for (int i = 0; i < n; i ++)
        ct[i] = c[i] - s[i];

    PDMatrix *pm2 = new PDMatrix(AT, n, m);
    pm2->Adjoint(cc, ct);
    SolverCG *so2 = new SolverCG(0.001, 0.01, 50, pm2);
    so2->CheckCgsOp();
    so2->CgsWy(v, v, cc);

    delete [] AT;
    delete [] ct;
    delete [] cc;
    delete [] bb;
}

int main (int argc, char **argv) {
   
    double start = time(NULL);
    //srand((unsigned)time(NULL));

    int n = 3;
    int m = 2;
    int ng = n * n;
    int na = m * n;

    double *A = new double[na];   // size(A) = m x n
    double *b = new double[m];
    double *c = new double[n];

    // get RHS of equation
    SetPara(A, b, c, na);


    // set a solution
    double *x = new double[n];
    double *v = new double[m];
    double *s = new double[n];
    // SetSolution(x, v, s, A, b, c, n, m);

        
    // define the specific problem
    KKTlp *kktm = new KKTlp(A, m, n);
    //kktm->PrintMatrix(A, m, n, "Aopr");
    //kktm->PrintVector(x, n, "xxx");

    // define the solver
    SolverLin *sollp = new SolverLin (kktm);
    

    //sollp->GetRhsBC(b, c, x, v, s);
    //kktm->PrintVector(b, m, "bbb");
    //kktm->PrintVector(c, n, "ccc");

    kktm->PrintMatrix(A, m, n, "Aopr");
    kktm->PrintVector(b, m, "bbb");
    kktm->PrintVector(c, n, "ccc");
	    

    // solve the problem
    double *x0 = new double[n];
    double *v0 = new double[m];
    double *s0 = new double[n];
    sollp->InnerPoint(x0, v0, s0, b, c);

    // print result
    cout << "index   predicted   true " << endl;
    for (int i = 0; i < n; i ++) {
        cout << " i, x = " << i << ", " << x0[i] << ", " << x[i] << endl; 
    }
    for (int i = 0; i < m; i ++) {
        cout << " i, v = " << i << ", " << v0[i] << ", " << v[i] << endl; 
    }
    for (int i = 0; i < n; i ++) {
        cout << " i, s = " << i << ", " << s0[i] << ", " << s[i] << endl; 
    }

    double stop = time(NULL);
    double durationTime = (double)difftime(stop, start);
    cout << endl << "Cost：" << durationTime << " sec" << endl;

    kktm->Free();
}
