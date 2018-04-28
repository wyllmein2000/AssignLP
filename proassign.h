#ifndef PROASSIGN_H
#define PROASSIGN_H
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
/*
#include <stdlib.h>
#include <set>
#include <map>
*/
#include "wyarray.h"
using namespace std;

//class KKTlp: public KKTEcqFunc {
//class KKTlp: public KKTMatrix {
class ProAssign {

    public:
	int nusr, nmsg;
	int nx, ne, ni;
        double lambda;  
	double *w, *bot;

	string regulation;

        ProAssign (double, double *, double *, int, int, string);

	void GetInitPoint (double *, double *, double *);
        void ComputeConstraint(double *, double *, double *);
        void ComputeSlash(double *, double *, double *, double);
	void ComputeGradient(double *, double *, double *, double *, double *, double *, double *, double);

        double Loss (double *);
        double LossAug (double *, double *, double *, double *, double *, double *, double);

	void round(double *, double *, int);
	void flow(double *, double *);
	double yield(double *);
	double misfit(double *);

        void PrintVector(double *, int, char *str);
        void PrintMatrix(double *, int, int, char *str);
        void WriteMatrix(double *, int, int, char *str);
	void printResult(string, double *, double *, int);
	void Free ();
};
#endif
