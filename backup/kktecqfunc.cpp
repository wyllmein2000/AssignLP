#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <set>
#include <map>
#include "wyarray.h"
#include "kktecqfunc.h"
using namespace std;


/* -------------------------------------------------------------------
 * KKT condition: 
 *   T(A) * x + s = c
 *      A * v     = w
 *   v[i] * s[i]  = 0
 *     (v, s)  >= 0
 *
 * KKT matrix:
 * | -Ds*Dv  Dv*T(A) |
 * |                 | dz = r
 * |   A*Dv    0     |
 * where dz = T(dv, dx)
 *
 * x = (v, x, s)
 * * dimension
 * size(A) = (nx, ns)
 * size(x) = nx
 * size(s) = ns
 * size(v) = nv
   -----------------------------------------------------------------*/
KKTEcqFunc::KKTEcqFunc (double *score, double *bottom, double *upper, int nusr, int nmsg) {

   this->nusr = nusr;
   this->nmsg = nmsg;
   this->nx = nusr * nmsg;
   this->ns = nusr + 2 * nmsg + nx;
   this->nv = this->ns;
   this->nm = this->nv + this->nx + this->ns;

   this->x0 = new double[nx];
   this->v0 = new double[nv];
   this->s0 = new double[ns];
   this->xo = new double[nm];

   this->pf = new double[nv];
   this->pa = new double[nv];

   this->w = new double[nx];
   memcpy(this->w, score, nx * sizeof(double));

   this->c = new double[ns];
   memset(this->c, 0, ns * sizeof(double));

   for (int i = 0; i < nusr; i ++)
       this->c[i] = 1.0;
   for (int i = 0; i < nmsg; i ++)
       this->c[nusr + i] = -bottom[i];
   for (int i = 0; i < nmsg; i ++)
       this->c[nusr + nmsg + i] = upper[i];

   this->init();
}

  
void KKTEcqFunc::init () {

   /*
   cout << " --- parameters --- " << endl;
   cout << " lambda = " << lambda << endl;
   cout << "  gamma = " <<  gamma << endl;
   */

   // initial value (v0, x0, s0)
   // (v0)
   for (int i = 0; i < nv; i ++) {
	v0[i] = 0.001;
   }
   // (x0)
   for (int i = 0; i < nusr; i ++) {
       for (int j = 0; j < nmsg; j ++) {
           x0[i * nmsg + j] = 1.0 / nmsg;
           //x0[i * nmsg + j] = 0.1 * (1 + rand() % 9);
       }
   }
   // (s0)
   for (int i = 0; i < ns; i ++) {
	s0[i] = 0.001;
   }
   
   // memcpy(xo, x0, nm * sizeof(double));
}


void KKTEcqFunc::free() {
   delete [] this->w;
   delete [] this->c;

   delete [] this->x0;
   delete [] this->v0;
   delete [] this->s0;
   delete [] this->xo;
   delete [] this->pf;
   delete [] this->pa;
}



void KKTEcqFunc::UpdateX (double *x) {
   memcpy(this->x0, x, nm * sizeof(double));
}

void KKTEcqFunc::GetX (double *x) {
   memcpy(x, this->v0, nv * sizeof(double));
   memcpy(&x[nv], this->x0, nx * sizeof(double));
   memcpy(&x[nv+nx], this->s0, ns * sizeof(double));
}

void KKTEcqFunc::GetDim (int *m) {
   *m = this->nm;
}

/* ----------------------------
 * | -DsDv   DvA* |   |dv|
 * |              | * |  | = r
 * |  A*Dv   0    |   |dx|
 * ---------------------------*/

/* Ax = y */
void KKTEcqFunc::Forward (double *y, double *x) {
    int n1 = nusr + nmsg;
    int n2 = nusr + nmsg * 2;
    for (int i = 0; i < this->nv; i ++)
	pf[i] = x[i] * v0[i];

    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[k + j] += pf[i] - pf[nusr + j] + pf[n1 + j] - pf[n2 + k + j];
	}
    }
}

/* y = A'x */
void KKTEcqFunc::Adjoint (double *y, double *x) {
    int n1 = nusr + nmsg;
    int n2 = nusr + nmsg * 2;
    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[i] += x[k + j];
	    y[nusr + j] -= x[k + j];
	    y[n1 + j] += x[k + j];
	    y[n2 + k + j] -= x[k + j];
	}
    }
    for (int i = 0; i < nv; i ++)
        y[i] *= v0[i];
}

void KKTEcqFunc::op(double *y, double *x) {
    memset(pa, 0, nv * sizeof(double));
    this->Adjoint(pa, x);
    this->Forward(y, pa);
}

void KKTEcqFunc::GetVo () {

    // CG parameters
    int cg_niter = 50;
    double cg_eps = 1e-06;
    double cg_wns = 0.01;

    double *dx = new double[nx];
    memset(dx, 0, nx * sizeof(double));

    // solving AA'v = Ar

    // define a PD function
    PDMtrFunc pdf(nusr, nmsg);

    // RHS
    pdf.Init(v0, w);

    // define a CG operator used in computing dx
    SolverCG solcg (cg_eps, cg_wns, cg_niter, this->nv);

    // check the operator
    // solcg.CheckCgsOp(pdf);

    // CG
    solcg.CgsWy (v0, dx, w, pdf);

    pdf.Free();
    delete [] dx;
}

void KKTEcqFunc::P_opr(double *y, double *r) {

    // CG parameters
    int cg_niter = 50;
    double cg_eps = 1e-06;
    double cg_wns = 0.01;

    double *av = new double[nx];
    double *b0 = new double[nv];
    double *dx = new double[nx];
    memset(dx, 0, nx * sizeof(double));

    // solving AA'v = Ar

    // define a PD function
    PDMtrFunc pdf(nusr, nmsg);

    // RHS
    pdf.Forward(b0, r);
    pdf.Init(v0, b0);

    // define a CG operator used in computing dx
    SolverCG solcg (cg_eps, cg_wns, cg_niter, this->nv);

    // check the operator
    // solcg.CheckCgsOp(pdf);

    // CG
    solcg.CgsWy (dx, dx, b0, pdf);

    // y = r - A'v
    pdf.Adjoint(av, dx);
    VectorSub(y, r, av, nx);

    pdf.Free();
    delete [] av;
    delete [] b0;
    delete [] dx;
}

void KKTEcqFunc::G_opr(double *y, double *x) {

    int n1 = nv + nx;
    memset(y, 0, this->nm * sizeof(double));

    for (int i = 0; i < this->nv; i ++) {
        y[i] = v0[i] * s0[i] * x[i];
    }
}

/* ----------------------------
 * | df/dx + (A*)v |
 * |               | = - residual
 * | b - Ax        | 
 * ---------------------------*/
void KKTEcqFunc::residual_aff(double *y, double *x) {

    double *y_fwd = new double[nx];
    double *y_adj = new double[ns];

    memset(y_fwd, 0, nx * sizeof(double));
    memset(y_adj, 0, ns * sizeof(double));

    this->Forward(y_fwd, x);
    this->Adjoint(y_adj, &x[ns]);

    for (int i = 0; i < ns; i ++) {
	y[i] = v0[i] * (y_adj[i] - c[i]);
    }

    for (int i = 0; i < nx; i ++) {
        y[ns + i] = w[i] - y_fwd[i];	
    }

    delete [] y_fwd;
    delete [] y_adj;
}

void KKTEcqFunc::maxres(double *r1, double *r2, double *x) {
    *r1 = vector_amax(x, nusr * nmsg);
    *r2 = vector_amax(&x[nusr * nmsg], nusr);
}

void KKTEcqFunc::round(double *y, double *x) {
    memset(y, 0, nusr * nmsg * sizeof(double));
    for (int iusr = 0; iusr < nusr; iusr ++) {
	int ks = iusr * nmsg;
        int kmsg = vector_max_index (&x[ks], nmsg);
        for (int imsg = 0; imsg < nmsg; imsg ++) {
            if (imsg == kmsg) 
	       y[ks + imsg] = 1.0;
	    else
	       y[ks + imsg] = 0.0;
	}
    }
}


void KKTEcqFunc::flow(double *y, double *x) {
    memset(y, 0, nmsg * sizeof(double));
    for (int iusr = 0; iusr < nusr; iusr ++)
    for (int imsg = 0; imsg < nmsg; imsg ++)
        y[imsg] += x[iusr * nmsg + imsg];
}

double KKTEcqFunc::yield(double *x) {
    int ns = nusr * nmsg;
    double y = 0.0;
    for (int i = 0; i < ns; i ++) {
	y += w[i] * x[i];
    }
    return y / ns;
}

double KKTEcqFunc::entropy(double *x, int n) {
    int ns = nusr * nmsg;
    double a = 0.0;
    double *y = new double[nmsg];
    this->flow(y, x);
    for (int i = 0; i < nmsg; i ++) {
	if (y[i] > 1e-30) 
	   a += y[i] * log(y[i]);
    }
    delete [] y;
    return a / nmsg;
}

double KKTEcqFunc::misfit(double *x) {
    int ns = nusr * nmsg;
    double y = 0.0;
    for (int i = 0; i < nx; i ++) {
	y += w[i] * x[i];
    }
    return y;
}


void KKTEcqFunc::printResult(string outputFileName) {

    double e0, y0;
    double entropy, yield;
    double entropy_r, yield_r;
    double r1, r2;

    double *xr = new double[nm];
    double *a0 = new double[nmsg];
    double *actual = new double[nmsg];
    double *actual_r = new double[nmsg];

    this->flow(a0, this->xo);
    e0 = this->entropy(a0, nmsg);
    y0 = this->yield(this->xo);

    this->round(xr, this->x0);
    this->flow(actual, this->x0);
    entropy = this->entropy(actual, nmsg);
    yield = this->yield(this->x0);

    this->flow(actual_r, xr);
    entropy_r = this->entropy(actual_r, nmsg);
    yield_r = this->yield(xr);

    cout << endl;
    cout << " --- final x ---- " << endl;
    //exaio.printArray(x1, nm);

    ofstream fp(outputFileName, ios::app);
   // if (fp) {
      fp << " nusr = " << nusr << endl;
      fp << " nmsg = " << nmsg << endl;
      fp << " lambda = " << lambda << endl;
      fp << " gamma = " << gamma << endl;
      fp << " flow for each msg: " << endl;
      for (int i = 0; i < nmsg; i ++) {
	  fp << setw(10) << i << setw(10) << a0[i] << " " << actual[i] << " " << actual_r[i] << " " << c[nusr+i] << " " << c[nusr+nmsg+i] << endl;
      }
      fp << setw(10) << "yield " << setw(10) << y0 << " " << yield << " " << yield_r << endl;
      fp << setw(10) << "entropy " << setw(10) << e0 << " " << entropy << " " << entropy_r << endl;
      fp << endl;
      fp.close();
  //  }

    delete [] xr;
    delete [] a0;
    delete [] actual;
    delete [] actual_r;
}
