#include <cstring>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "assign.h"
using namespace std;

assign::assign (double lambda, double gamma, double *score, double *target, double *ass, int nusr, int nmsg) {
   this->lambda = lambda;
   this->gamma = gamma;
   this->target = target;
   this->score = score;
   this->ass = ass;
   this->nusr = nusr;
   this->nmsg = nmsg;
}

/* ----------------------------
 * | H   A* |
 * |        |x = y
 * | A   0  | 
 * ---------------------------*/

/* Ax = y */
void assign::forward(double *y, double *x) {
    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[i] += x[k + j];
	}
    }
}

/* (A*)x = y */
void assign::adjoint(double *y, double *x) {
    for (int i = 0; i < nusr; i ++) {
	int k = i * nmsg;
	for (int j = 0; j < nmsg; j ++) {
	    y[k + j] += x[i];
	}
    }

}

void assign::op(double *y, double *x) {
    int n = nusr * nmsg;
    int m = nusr;

    memset(y, 0, (n + m) * sizeof(double));

    for (int iusr = 0; iusr < nusr; iusr ++)
    for (int imsg = 0; imsg < nmsg; imsg ++) {
        int is = iusr * nmsg + imsg;
	y[is] += x[is] * lambda / ass[is];
	
	for (int jusr = 0; jusr < nusr; jusr ++) {
	    int js = jusr * nmsg + imsg;
	    y[is] += gamma * x[js];
	}
    }

    this->forward(&y[n], x);
    this->adjoint(y, &x[n]);
}

/* ----------------------------
 * | df/dx + (A*)v |
 * |               | = - residual
 * | b - Ax        | 
 * ---------------------------*/
void assign::residual(double *y, double *x) {

    int ns = nusr * nmsg;
    double *x_sum = new double[nmsg];
    double *y_fwd = new double[nusr];
    double *y_adj = new double[ns];

    memset(x_sum, 0, nmsg * sizeof(double));
    memset(y_fwd, 0, nusr * sizeof(double));
    memset(y_adj, 0, ns   * sizeof(double));

    this->forward(y_fwd, x);
    this->adjoint(y_adj, &x[ns]);

    for (int iusr = 0; iusr < nusr; iusr ++) {
         int ks = iusr * nmsg;
         for (int imsg = 0; imsg < nmsg; imsg ++) {
              x_sum[imsg] += x[ks + imsg];
	 }
    }

    for (int iusr = 0; iusr < nusr; iusr ++) {

        int ks = iusr * nmsg;

        for (int imsg = 0; imsg < nmsg; imsg ++) {
             int k = ks + imsg;
	     y[k] = score[k] - lambda * (1.0 + log(x[k])) - gamma * (x_sum[imsg] - target[imsg]) - y_adj[k];
	        //y[k] = score[k] - gamma * (x_sum[imsg] - target[imsg]) - y_adj[k];
	}

	    // cout << " iusr=" << iusr << ", yfwd=" << y_fwd[iusr] << endl;
	y[ns + iusr] = 1.0 - y_fwd[iusr];
    }

    delete [] x_sum;
    delete [] y_fwd;
    delete [] y_adj;
}

void assign::flow(double *y, double *x) {
    memset(y, 0, nmsg * sizeof(double));
    for (int iusr = 0; iusr < nusr; iusr ++)
    for (int imsg = 0; imsg < nmsg; imsg ++)
        y[imsg] += x[iusr * nmsg + imsg];
}

double assign::yield(double *x) {
    int ns = nusr * nmsg;
    double y = 0.0;
    for (int i = 0; i < ns; i ++) {
	y += score[i] * x[i];
    }
    return y;
}

double assign::misfit(double *x) {
    int ns = nusr * nmsg;
    double y0 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    for (int i = 0; i < ns; i ++) {
	y0 += score[i] * x[i];
        y1 += x[i] * log(x[i]);
    }
    for (int imsg = 0; imsg < nmsg; imsg ++) {
	double yt = -target[imsg];
	for (int iusr = 0; iusr < nusr; iusr ++) {
	    yt += x[iusr * nmsg + imsg];
	}
        y2 += yt * yt;
    }
    return y0 - lambda * y1 - 0.5 * gamma * y2;
}
