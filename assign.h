#ifndef ASSIGN_H
#define ASSIGN_H

class assign {
    private:
        double lambda;
	double gamma;
	double *score;
	double *target;
	double *ass;
	int nusr;
	int nmsg;

    public:
        assign (double lambda, double gamma, double *score, double *target, double *ass, int nusr, int nmsg);
        void forward(double *y, double *x);
        void adjoint(double *y, double *x);
        void op(double *y, double *x);
        void residual(double *y, double *x);
        void flow(double *y, double *x);
        double yield(double *x);
        double misfit(double *x);
}
;


#endif
