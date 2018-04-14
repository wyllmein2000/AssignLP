#include "proassign.h"
using namespace std;


/* -------------------------------------------------------------------
   -----------------------------------------------------------------*/
ProAssign::ProAssign (double *score, double *bottom, int nusr, int nmsg) {

   this->lambda = 0.0;

   this->nusr = nusr;
   this->nmsg = nmsg;
   this->nx = nusr * nmsg;
   this->ne = nusr;
   this->ni = nmsg;

   this->w = new double[nx];
   memcpy(this->w, score, nx * sizeof(double));

   this->bot = new double[ni];
   memcpy(this->bot, bottom, ni * sizeof(double));

}

  
void ProAssign::GetInitPoint (double *x, double *ce, double *cn) {

   for (int i = 0; i < nx; i ++)
       x[i] = 1.0 / this->nmsg;

   /*
   for (int i = 0; i < nusr; i++) {
       int k = i * nmsg;
       for (int j = 0; j < nmsg; j ++) {
	   x[k + j] = (1.66 * bot[j]) / nusr;
       }
   }*/
   for (int i = 0; i < ne; i ++)
       ce[i] = 0.1;	
   for (int j = 0; j < ni; j ++)
       cn[j] = 0.1;	
}



/* ve --- equality constraint 
 * vn --- inequaltiy constraint 
 */
void ProAssign::ComputeConstraint(double *ve, double *vi, double *x) {
   int k;
   for (int i = 0; i < nusr; i ++)
       ve[i] = -1.0;

   for (int j = 0; j < nmsg; j ++)
       vi[j] = -1.0 * this->bot[j];

   for (int i = 0; i < nusr; i ++) {
       k = i * nmsg;
       for (int j = 0; j < nmsg; j ++) {
	   ve[i] += x[k + j];
	   vi[j] += x[k + j];
       }
   }
}

void ProAssign::ComputeSlash(double *s, double *cn,  double *vn, double sigma) {
   for (int j = 0; j < nmsg; j ++) {
       s[j] = vn[j] - cn[j] / sigma;
       if (s[j] < 0) s[j] = 0.0;
   }
}

void ProAssign::ComputeGradient(double *g, double *x, double *ce, double *cn, double *ve, double *vn, double *s, double sigma) {
   int k, m;
   for (int i = 0; i < nusr; i ++) {
       k = i * nmsg;
       for (int j = 0; j < nmsg; j ++) {
	   m = k + j;
           g[m] = -w[m] - ce[i] - cn[j] + sigma * (ve[i] + vn[j] - s[j]);
	   // maximum entropy
	   // g[m] += this->lambda * (1.0 + log(x[m]));
	   // minimum energy
	   g[m] += this->lambda * x[m];
       }
   }
}

/*
void ProAssign::UpdateX(double *x, double *dx, double alpha) {
   for (int i = 0; i < this->nx; i++) {
       x[i] -= alpha * dx[i];
       if (x[i] < 0.0) x[i] = 0.0;
   }
}*/

double ProAssign::Loss(double *x) {
   double y = 0.0;
   for (int i = 0; i < this->nx; i ++)
       y += -1.0 * this->w[i] * x[i] + 0.5 * this->lambda * x[i] * x[i];
   return y;
}

double ProAssign::LossAug(double *x, double *ce, double *cn, double *ve, double *vn, double *s, double sigma) {
   double y = this->Loss(x);
   double z = 0.0;
   for (int i = 0; i < this->ne; i ++) {
       y -= ce[i] * ve[i];
       z += ve[i] * ve[i];
   }
   for (int i = 0; i < this->ni; i ++) {
       y -= cn[i] * (vn[i] - s[i]);
       z += (vn[i] - s[i]) * (vn[i] - s[i]);
   }
   y += 0.5 * sigma * z;
   return y;
}



void ProAssign::round(double *y, double *x) {
    memset(y, 0, nusr * nmsg * sizeof(double));
    srand((unsigned)time(NULL));
    this->WriteMatrix(x, nusr, nmsg, "output/msg_open_log");
    int flag = 1;
    if (flag == 0) {
    for (int iusr = 0; iusr < nusr; iusr ++) {
	int ks = iusr * nmsg;
        //int kmsg = vector_max_index (&x[ks], nmsg);
        int kmsg = vector_max_index (&y[ks], nmsg);
        for (int imsg = 0; imsg < nmsg; imsg ++) {
            if (imsg == kmsg) 
	       y[ks + imsg] = 1.0;
	    else
	       y[ks + imsg] = 0.0;
	}
    }
    }
    else {
    for (int iusr = 0; iusr < nusr; iusr ++) {
	int ks = iusr * nmsg;
	int imsg = 0;
	while (ks >= 0) {
	    float r = rand()/(RAND_MAX + 1.0);
	    if (imsg >= nmsg) 
	       imsg = imsg - nmsg;
	    if (r < x[ks + imsg]) {
               y[ks + imsg] = 1.0;
	       break;
	    }
	    imsg += 1;
	}
    }
    }

    this->WriteMatrix(y, nusr, nmsg, "output/msg_open_log_r");
}


void ProAssign::flow(double *y, double *x) {
    memset(y, 0, nmsg * sizeof(double));
    for (int iusr = 0; iusr < nusr; iusr ++)
    for (int imsg = 0; imsg < nmsg; imsg ++)
        y[imsg] += x[iusr * nmsg + imsg];
}

double ProAssign::yield(double *x) {
    double y = 0.0;
    for (int i = 0; i < this->nx; i ++) {
	y += this->w[i] * x[i];
    }
    return y / this->nx;
}

/*
double ProAssign::entropy(double *x, int n) {
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
}*/


void ProAssign::PrintVector(double *x, int m, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) 
       cout << " i = " << i << ", " << x[i] << endl;
}
  
void ProAssign::PrintMatrix(double *x, int m, int n, char *str) {
   cout << endl << str << endl;
   for (int i = 0; i < m; i ++) {
       cout << " i=" << i << " "; 
       for (int j = 0; j < n; j ++)
           cout << x[i*n+j] << " ";
       cout << endl;
   }
}
  
void ProAssign::WriteMatrix(double *x, int m, int n, char *str) {
   ofstream fp(str);
   for (int i = 0; i < m; i ++) {
       fp << i << " "; 
       for (int j = 0; j < n; j ++)
           fp << x[i*n+j] << " ";
       fp << endl;
   }
   fp.close();
}
  
void ProAssign::printResult(string outputFileName, double *x0, double *x) {

    double e0, y0;
    double entropy, yield;
    double entropy_r, yield_r;
    double s_x0, s_x, s_xr;
    double r1, r2;

    double *xr = new double[nusr * nmsg];
    double *a0 = new double[nmsg];
    double *actual = new double[nmsg];
    double *actual_r = new double[nmsg];
    double *avsc = new double[nmsg];

    // init
    this->flow(a0, x0);
    //e0 = this->entropy(a0, nmsg);
    y0 = this->yield(x0);
    s_x0 = vector_sum(a0, nmsg);

    // solved x
    this->flow(actual, x);
    //entropy = this->entropy(actual, nmsg);
    yield = this->yield(x);
    s_x = vector_sum(actual, nmsg);


    // round
    this->round(xr, x);
    this->flow(actual_r, xr);
    //entropy_r = this->entropy(actual_r, nmsg);
    yield_r = this->yield(xr);
    s_xr = vector_sum(actual_r, nmsg);

    // average score
    memset(avsc, 0, nmsg * sizeof(double));
    for (int i = 0; i < nusr; i++) 
    for (int j = 0; j < nmsg; j++) 
        avsc[j] += this->w[i * nmsg + j];

    cout << endl;
    cout << " --- final x ---- " << endl;
    //exaio.printArray(x1, nm);

    //ofstream fp(outputFileName.c_str(), ios::app);
    ofstream fp(outputFileName, ios::app);
   // if (fp) {
      fp << " nusr = " << nusr << endl;
      fp << " nmsg = " << nmsg << endl;
      fp << " lambda = " << lambda << endl;
      fp << " flow for each msg: " << endl;
      for (int i = 0; i < nmsg; i ++) {
	  fp << setw(10) << i << setw(10) << a0[i] << setw(10) << actual[i] << setw(10) << actual_r[i] << setw(10) << this->bot[i] << setw(10) << avsc[i] << endl;
      }
      fp << setw(10) << "sum" << setw(10) << s_x0 << setw(10) << s_x << setw(10) << s_xr << setw(10) << " " << setw(10) << " " << endl;
      fp << setw(10) << "yield " << setw(10) << y0 << " " << yield << " " << yield_r << endl;
      //fp << setw(10) << "entropy " << setw(10) << e0 << " " << entropy << " " << entropy_r << endl;
      fp << endl;
      fp.close();
  //  }

    delete [] xr;
    delete [] a0;
    delete [] actual;
    delete [] actual_r;
}


void ProAssign::Free() {
   delete [] this->w;
   delete [] this->bot;
}

