#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <map>
#include "inout.h"
using namespace std;

inout::inout () {
}

void inout::readtext_example (double *score, int nuser, int nmsg, int minMsgId) {
   ifstream inf;
   inf.open("/Users/wuyan/workdir/aaa.txt", ifstream::in);

   const int cnt = 4;
   size_t icomma1 = 0;
   size_t icomma2 = 0;
   string str, user1, user2;
   int msgId, iusr = 0, k = 0;
   double point;

   string line;

   while (!inf.eof()) {
      getline(inf, line);

      icomma1 = line.find(',', 0);
      user1 = line.substr(0, icomma1);
      //cout << "j=0:" << user1.c_str() << endl;

      for (int j = 1; j < cnt; j ++) {
	  icomma2 = line.find(',', icomma1 + 1);
          str = line.substr(icomma1 + 1, icomma2 - icomma1 - 1);
          //cout << "j=" << j << ":" << str.c_str() << endl;
	  if (j == 1)
             msgId = atoi(str.c_str()) - minMsgId;
	  else if (j == 2)
             point = atof(str.c_str());	  
	  icomma1 = icomma2;
      }

      if (k > 1 && strcmp(user1.c_str(), user2.c_str()) != 0) iusr += 1;

      int index = iusr * nmsg + msgId;
      if (index >= nuser * nmsg || index < 0) {
	 cout << " Error: out of range ..." << endl;
	 break;
      }
      else {
         score[index] = point;
      }

      user2 = user1;
      k += 1;
   }


}






void inout::printArray(double *x, int n) {
   for (int i = 0; i < n; i ++) {
       cout << "[" << i << "]=" << x[i] << endl;
   }
}

void inout::print2Array(double *x, double *y, int n) {
   for (int i = 0; i < n; i ++) {
       cout << "i = " << i << ", x=" << x[i] << ", y=" << y[i] << endl;
   }
}

void inout::printArray2d(double *x, int n, int m) {
   for (int i = 0; i < n; i ++) {
       cout << endl;
       for (int j = 0; j < m; j ++) {
	   cout << "[" << i << "][" << j << "]=" << x[i*m+j] << endl;
       }
   }
}


void inout::writedata (char *fn, double *p, int n) {
     FILE *fp = fopen(fn,"wb");
     for (int i=0; i<n; i++)
         fwrite(&p[i],sizeof(double),1,fp);
     fclose(fp);
}
