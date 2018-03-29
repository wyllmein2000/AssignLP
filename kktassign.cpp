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
#include "kktassign.h"
using namespace std;


/* -------------------------------------------------------------------
 * KKT condition: 
 *   T(A) * x + s = c
 *      A * v     = w
 *   v[i] * s[i]  = 0
 *     (v, s)  >= 0
 *
 *
 * | Ds*Dv  Dv*A' |   |-dv|   |Rxs - Dv*Rc|
 * |              | * |   | = |           |
 * |  A*Dv    0   |   | dx|   |    Rb     |
 *
 * where Rb = Av - w 
 *       Rc = A'x + s - c
 *       Rxs = Dv*Ds*e
 *
 *
 * KKT matrix:
 * | Ds*Dv  Dv*A' |   |-dv|
 * |              | * |   | = r
 * |  A*Dv    0   |   | dx|
 * and ds = -Ds*(e+dv)
 *
 * x = (v, x, s)
 * * dimension
 * size(A) = (nx, ns)
 * size(x) = nx
 * size(s) = ns
 * size(v) = nv
   -----------------------------------------------------------------*/
KKTAssign::KKTAssign (double *score, double *bottom, double *upper, int nusr, int nmsg) {

   this->nusr = nusr;
   this->nmsg = nmsg;
   this->nx = nusr * nmsg;
   this->ns = nusr + 2 * nmsg + nx;
   this->nv = this->ns;
   this->nm = this->nv + this->nx;
   this->m = ns;
   this->n = nx;

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

  
void KKTAssign::init () {

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


void KKTAssign::free() {
   delete [] this->w;
   delete [] this->c;

   delete [] this->x0;
   delete [] this->v0;
   delete [] this->s0;
   delete [] this->xo;
   delete [] this->pf;
   delete [] this->pa;
}



void KKTAssign::UpdateX (double *x) {
   memcpy(this->x0, x, nm * sizeof(double));
}

void KKTAssign::GetX (double *x) {
   memcpy(x, this->v0, nv * sizeof(double));
   memcpy(&x[nv], this->x0, nx * sizeof(double));
   memcpy(&x[nv+nx], this->s0, ns * sizeof(double));
}

void KKTAssign::GetDim (int *m) {
   *m = this->nm;
}

/* ds = -Ds * (e + dv) */
void KKTAssign::ComputeDs(double *s, double *v) {
   for (int i = 0; i < m; i++)
       s[i] = -s0[i] * (1.0 + v[i]);
}


void KKTAssign::GetRhsPredictor (double *ru, double *rb) {
   double *rc = new double[m];
   this->Adjoint(rc, x0);
   for (int i = 0; i < m; i++) {
       rc[i] = rc[i] + s0[i] - c[i];
       ru[i] = v0[i] * s0[i] - v0[i] * rc[i];
   }
   this->Forward(rb, v0);
   for (int i = 0; i < n; i++) {
       rb[i] = w[i] - rb[i];
   }
}


void KKTAssign::GetRhsCorrector (double *ru, double *rb, double *dr) {
   double *rc = new double[m];
   this->Adjoint(rc, x0);
   for (int i = 0; i < m; i++) {
       rc[i] = rc[i] + s0[i] - c[i];
       ru[i] = v0[i] * s0[i] - v0[i] * rc[i] + dr[i];
   }
   this->Forward(rb, v0);
   for (int i = 0; i < n; i++) {
       rb[i] = w[i] - rb[i];
   }
}



/* Ax = y */
void KKTAssign::Forward (double *y, double *x) {
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
void KKTAssign::Adjoint (double *y, double *x) {
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


void KKTAssign::G_opr(double *y, double *x) {
    memset(y, 0, this->nv * sizeof(double));
    for (int i = 0; i < this->nv; i ++) {
        y[i] = v0[i] * s0[i] * x[i];
    }
}


void KKTAssign::Opr(double *y, double *x, int flag) {
    if (flag == 1) {
       double *p = new double[m];
       memset(p, 0, m * sizeof(double));
       memset(y, 0, n * sizeof(double));
       this->Forward(p, x);
       this->Adjoint(y, p);
       delete [] p;
    }
    else if (flag == 0) {
       double *p = new double[n];
       memset(p, 0, n * sizeof(double));
       memset(y, 0, m * sizeof(double));
       this->Adjoint(p, x);
       this->Forward(y, p);
       delete [] p;
    }
}



void KKTAssign::KKTfunc(double *yc, double *yb, double *x, double *v) {
    double *p = new double[n];
    memset(p, 0, n * sizeof(double));
    memset(yb, 0, m * sizeof(double));
    this->Forward(yb, x);
    this->Adjoint(p, v);
    this->G_opr(yc, x);
    VectorSub(yc, p, yc, n);
    delete [] p;
}



















/* ----------------------------
 * | df/dx + (A*)v |
 * |               | = - residual
 * | b - Ax        | 
 * ---------------------------*/
void KKTAssign::residual_aff(double *y, double *x) {

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

void KKTAssign::maxres(double *r1, double *r2, double *x) {
    *r1 = vector_amax(x, nusr * nmsg);
    *r2 = vector_amax(&x[nusr * nmsg], nusr);
}

void KKTAssign::round(double *y, double *x) {
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


void KKTAssign::flow(double *y, double *x) {
    memset(y, 0, nmsg * sizeof(double));
    for (int iusr = 0; iusr < nusr; iusr ++)
    for (int imsg = 0; imsg < nmsg; imsg ++)
        y[imsg] += x[iusr * nmsg + imsg];
}

double KKTAssign::yield(double *x) {
    int ns = nusr * nmsg;
    double y = 0.0;
    for (int i = 0; i < ns; i ++) {
	y += w[i] * x[i];
    }
    return y / ns;
}

double KKTAssign::entropy(double *x, int n) {
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

double KKTAssign::misfit(double *x) {
    int ns = nusr * nmsg;
    double y = 0.0;
    for (int i = 0; i < nx; i ++) {
	y += w[i] * x[i];
    }
    return y;
}


void KKTAssign::printResult(string outputFileName) {

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
