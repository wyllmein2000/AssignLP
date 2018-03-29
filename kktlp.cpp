#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include "wyarray.h"
#include "kktlp.h"
using namespace std;


/* -------------------------------------------------------------------
 * min(c'x)  s.t. Ax = b, x >= 0
 * max(b'v)  s.t. A'v + s = c, s >= 0
 *
 * KKT condition
 * A'v + s = c
 * Ax = b
 * x[i] * s[i] = 0
 * (x, s) >= 0
 *
 * KKT matrix:
 * | 0  A' I|   |dx|    |Rc |
 * | A  0  0| * |dv| = -|Rb |
 * | Ds 0 Dx|   |ds|    |Rxs|
 * where Rc = A'v + s -c
 *       Rb = Ax - b
 *       Rxs = Dx*Ds*e
 *
 * In the corrector step
 * Rxs = (Dx * Ds + D(dx_aff) * D(ds_aff) - sigma * mu)*e
 *
 *
 * Reduced KKT matrix:
 * | DsDx  Dx*A'|   |-dz|    |Rxs-Dx*Rc|   | resu|
 * |            | * |   | =  |         | = |     |
 * | A*Dx   0   |   | dv|    |-(-Rb)   |   |-resb|
 * and ds = -Ds(e+dz)
 *     resu = Rxs - Dx * Rc = Dx(c-A'v)
 *     resb = -Rb = b-Ax
 *
 * where G=Ds*Dx may be indefinite
 *       dx = Dx * dz
 *
 *
 * size(A) = (m, n)
 * size(G) = (n, n)
 * size(x) = size(s) = size(c) = n
 * size(v) = size(b) = m
   -----------------------------------------------------------------*/
KKTlp::KKTlp (double *a, int m, int n) {

   this->m = m;
   this->n = n;
   this->na = m * n;

   this->a = new double[na];
   memcpy(this->a, a, na * sizeof(double));

   this->x0 = new double[n];
   this->v0 = new double[m];
   this->s0 = new double[n];

   for (int i = 0; i < n; i ++) {
       this->x0[i] = 1.0;
       this->s0[i] = 1.0;
   }
   for (int i = 0; i < m; i ++) {
       this->v0[i] = 1.0;
   }

}


/* y = A D(x) x */
void KKTlp::Forward (double *y, double *x) {
    double *p = new double[n];
    for (int i = 0; i < n; i ++)
        p[i] = this->x0[i] * x[i];

    for (int i = 0; i < m; i ++) {
	int k = i * n;
	for (int j = 0; j < n; j ++) {
	    y[i] += a[k + j] * p[j];
	}
    }
    delete [] p;
}

/* y = D(x) A' x */
void KKTlp::Adjoint (double *y, double *x) {
    for (int i = 0; i < m; i ++) {
	int k = i * n;
	for (int j = 0; j < n; j ++) {
	    y[j] += a[k + j] * x[i];
	}
    }
    for (int i = 0; i < m; i ++)
        y[i] *= this->x0[i];
}

/* y = G x = D(s)D(x) x */
void KKTlp::G_opr(double *y, double *x) {
    memset(y, 0, n * sizeof(double));
    for (int i = 0; i < n; i ++) {
        y[i] = this->x0[i] * this->s0[i] * x[i];
    }
}


/* flag = 1,  dz = P * dx         *
 * flag = 0,  dx = inv(P) * dz    */
/* Here P = inv(x0)               */
void KKTlp::Precon(double *dx, double *dv, double *ds, int flag) {
    if (flag == 1) {
       for (int i = 0; i < n; i ++) {
	   dx[i] = dx[i] / this->x0[i];
       }
    }
    else if (flag == 0) {
       for (int i = 0; i < n; i ++) {
	   dx[i] = dx[i] * this->x0[i];
       }
    }
}




/* ifwd = 1, y = D(x)' A' A D(x) x */
/* ifwd = 0, y = A D(x) D(x)' A' x */
void KKTlp::Opr(double *y, double *x, int flag) {
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





/* yc = A'v - Gx 
 * yb = Ax        */
void KKTlp::KKTfunc(double *yc, double *yb, double *x, double *v) {
    double *p = new double[n];
    memset(p, 0, n * sizeof(double));
    memset(yb, 0, m * sizeof(double));
    this->Forward(yb, x);
    this->Adjoint(p, v);
    this->G_opr(yc, x);
    VectorSub(yc, p, yc, n);
    delete [] p;
}

void KKTlp::PrintVector(double *x, int m, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) 
       cout << " i = " << i << ", " << x[i] << endl;
}
  
void KKTlp::PrintMatrix(double *x, int m, int n, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) {
       cout << " i=" << i << " "; 
       for (int j = 0; j < n; j ++)
           cout << x[i*n+j] << " ";
       cout << endl;
   }
}
  


void KKTlp::UpdateX (double *x0, double *v0, double *s0) {
   memcpy(this->x0, x0, n * sizeof(double));
   memcpy(this->v0, v0, m * sizeof(double));
   memcpy(this->s0, s0, n * sizeof(double));
}
  

void KKTlp::Free() {
   delete [] this->a;
   delete [] this->x0;
   delete [] this->v0;
   delete [] this->s0;
}

