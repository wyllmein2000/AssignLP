#ifndef ASSIGN_H
#define ASSIGN_H
#include <string>
#include <iostream>
using namespace std;

class assign {
    private:
	int nusr;
	int nmsg;
	int ns;

        double lambda;
	double gamma;

	double *score;
	double *target;
	double *x0;
	double *xo;

	string fileName;

    public:

	int nm;

        assign (string fileName, int n);
        void readtext_usr_msg_score (string, int);

	void init ();
	void free ();
        void UpdateX (double *x);
        void GetX (double *x);
	void Constraint (double *x);

        void forward(double *y, double *x);
        void adjoint(double *y, double *x);
        void op(double *y, double *x);

        void residual(double *y, double *x);
        void maxres(double *r1, double *r2, double *x);

        void round(double *y, double *x);
        void flow(double *y, double *x);
        double entropy(double *x, int n);
        double yield(double *x);
        double misfit(double *x);

        void printResult(string outputFileName);
}
;


#endif
