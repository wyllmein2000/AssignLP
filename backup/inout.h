#ifndef INOUT_H
#define INOUT_H
#include <string>
using namespace std;

class inout {
   private:
      
   public:
      inout ();
      void readtext_example(double *, int, int, int);

      void printArray(double *x, int n);
      void print2Array(double *x, double *y, int n);
      void printArray2d(double *x, int n, int m);

      void writedata (char *fn, double *p, int n);

      void writeDataText();
      void printResult(double y, double ry, double e, double re, double *flow, double *rflow, double *target, int n);
};

#endif
