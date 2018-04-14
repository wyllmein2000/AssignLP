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
#include "assign.h"
using namespace std;

assign::assign (string fileName, int maxLine) {

   
   this->readtext_usr_msg_score (fileName, maxLine);

   //this->nusr = nusr;
   //this->nmsg = nmsg;
   this->ns = nusr * nmsg;
   this->nm = nusr + ns;

   this->target = new double[nmsg];
   this->x0 = new double[nm];
   this->xo = new double[nm];
}

  
void assign::init () {

   this->lambda = 1e-03;
   this->gamma = 1e-00;

   cout << " --- parameters --- " << endl;
   cout << " lambda = " << lambda << endl;
   cout << "  gamma = " <<  gamma << endl;

   //memset(score, 0, ns * sizeof(double));
   memset(target, 0, nmsg * sizeof(double));


   // read target ---  target flow for each msg
   // srand(seed);
   for (int j = 0; j < nmsg; j ++) {
	//target[j] = 60.0 + (rand() % 101 - 50);
	//target[j] = 3.0;
	target[j] = 50.0;
   }
   //target[0] = 4; target[1] = 4.0;

   // initial value (assignment)
   for (int i = 0; i < nusr; i ++) {
       for (int j = 0; j < nmsg; j ++) {
           x0[i * nmsg + j] = 1.0 / nmsg;
           //x0[i * nmsg + j] = 0.1 * (1 + rand() % 9);
       }
   }
   for (int i = 0; i < nusr; i ++) {
	x0[ns + i] = 0.001;
   }
   memcpy(xo, x0, nm * sizeof(double))

}

void assign::free() {
   delete [] this->score;
   delete [] this->target;
   delete [] this->x0;
}



// read data score & ????
void assign::readtext_usr_msg_score (string fileName, int maxLine) {
   ifstream inf;
   inf.open(fileName.c_str(), ifstream::in);
   cout << " read " << fileName.c_str() << endl;

   size_t index_separator;
   size_t index_separator_1;
   size_t index_separator_2;
   int iusr, imsg, usrNum, msgNum;
   int k = 0;
   double point;


   // read all usrId and msgId and save them in a set.
   set<string> usrIds;
   set<string> msgIds;
   string usrId, msgId;
   string str, line;


   while (!inf.eof()) {
      k += 1;
      if (k > maxLine) break;

      getline(inf, line);
      //if (k % 1000 == 0)
      //cout << line.c_str() << endl;
      //  cout << " k = " << k << endl;

      index_separator = line.find(',', 0);
      usrId = line.substr(0, index_separator);

      index_separator_1 = line.find(',', index_separator + 1);
      msgId = line.substr(index_separator + 1, index_separator_1 - index_separator - 1);

      if (!usrId.empty()) usrIds.insert(usrId);
      if (!msgId.empty()) msgIds.insert(msgId);

      // cout << usrId.c_str() << " " << msgId.c_str() << endl;
   }
   inf.close();

   usrNum = usrIds.size();
   msgNum = msgIds.size();

   cout << " number of users: " << usrNum << endl;
   cout << " number of messages: " << msgNum << endl;

   this->nusr = usrNum;
   this->nmsg = msgNum;

   map<int,string> id2usr, id2msg;
   map<string,int> usr2id, msg2id;
   k = 0;
   for (set<string>::iterator it = usrIds.begin(); it != usrIds.end(); it++) {
       id2usr.insert(pair<int,string>(k, *it));
       usr2id.insert(pair<string,int>(*it, k));
       k += 1;
   }
   k = 0;
   for (set<string>::iterator it = msgIds.begin(); it != msgIds.end(); it++) {
       id2msg.insert(pair<int,string>(k, *it));
       msg2id.insert(pair<string,int>(*it, k));
       k += 1;
   }

   this->score = new double[nusr * nmsg];

   k = 0;
   inf.open(fileName.c_str(), ifstream::in);

   while (!inf.eof()) {
      k += 1;
      if (k > maxLine) break;

      getline(inf, line);

      index_separator = line.find(',', 0);
      usrId = line.substr(0, index_separator);

      index_separator_1 = line.find(',', index_separator + 1);
      msgId = line.substr(index_separator + 1, index_separator_1 - index_separator - 1);

      index_separator_2 = line.find(',', index_separator_1 + 1);
      str = line.substr(index_separator_1 + 1, index_separator_2 - index_separator_1 - 1);

      if (!usrId.empty() && !msgId.empty()) {
	  iusr = usr2id[usrId];
          imsg = msg2id[msgId];
          int ks = iusr * nmsg + imsg;
	  if (ks < 0 || ks >= nusr * nmsg) continue;
          score[ks] = atof(str.c_str());
	  //cout << " ks = " << ks << ", score = " << score[ks] << endl;
      }

      // cout << usrId.c_str() << " " << msgId.c_str() << endl;
   }
   inf.close();

}

void assign::UpdateX (double *x) {
   memcpy(this->x0, x, nm * sizeof(double));
}

void assign::GetX (double *x) {
   memcpy(x, this->x0, nm * sizeof(double));
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
	y[is] += x[is] * lambda / x0[is];
	
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
    this->flow(x_sum, x);

    for (int iusr = 0; iusr < nusr; iusr ++) {

        int ks = iusr * nmsg;
        for (int imsg = 0; imsg < nmsg; imsg ++) {
             int k = ks + imsg;
	     y[k] = score[k] - lambda * (1.0 + log(x[k])) - gamma * (x_sum[imsg] - target[imsg]) - y_adj[k];
	}

	    // cout << " iusr=" << iusr << ", yfwd=" << y_fwd[iusr] << endl;
	y[ns + iusr] = 1.0 - y_fwd[iusr];
    }

    delete [] x_sum;
    delete [] y_fwd;
    delete [] y_adj;
}

void assign::maxres(double *r1, double *r2, double *x) {
    *r1 = vector_amax(x, nusr * nmsg);
    *r2 = vector_amax(&x[nusr * nmsg], nusr);
}

void assign::round(double *y, double *x) {
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

double assign::entropy(double *x, int n) {
    double y = 0.0;
    for (int i = 0; i < n; i ++) {
	if (x[i] > 1e-30) 
	   y += x[i] * log(x[i]);
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


void assign::printResult(string outputFileName) {

    double entropy, yield;
    double entropy_r, yield_r;
    double r1, r2;

    double *xr = new double[nm];
    double *actual = new double[nmsg];
    double *actual_r = new double[nmsg];

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
      fp << " nmsg = " << nusr << endl;
      fp << " lambda = " << lambda << endl;
      fp << " gamma = " << gamma << endl;
      fp << " flow for each msg: " << endl;
      for (int i = 0; i < nmsg; i ++) {
	  fp << setw(10) << i << setw(10) << actual[i] << " " << actual_r[i] << " " << target[i] << endl;
      }
      fp << setw(10) << "yield " << setw(10) << yield << " " << yield_r << endl;
      fp << setw(10) << "entropy " << setw(10) << entropy << " " << entropy_r << endl;
      fp << endl;
      fp.close();
  //  }

    delete [] xr;
    delete [] actual;
    delete [] actual_r;
}
