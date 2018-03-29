#include <fstream>
#include "solverlin.h"
/* -------------------------------------------------
 * algorithm
 *
 * ---------------------------------------------- */




/* solve linear programming with inner-point method */
/* -------------------------------------------------
 * min(c'x)  s.t. Ax = b, x >= 0
 * max(b'v)  s.t. A'v + s = c, s >= 0
 *
 * KKT condition
 * A'v + s = c
 * Ax = b
 * x[i] * s[i] = 0
 * (x, s) >= 0
 * ------------------------------------------------ */
SolverLin::SolverLin (KKTlp *ass) {
    this->inner = 10;
    this->outer = 5;
    this->niter = 5;
    this->mu = 0.6;
    this->eta = 0.9;

    this->iprint = 1;
    this->ass = ass;

    this->m = ass->m;
    this->n = ass->n;
}


void SolverLin::GetInitPoint (double *x, double *v, double *s, double *b, double *c, SolverPCG *solpcg) {

    double *p = new double[m];
    memset(p, 0, m * sizeof(double));
    memset(x, 0, n * sizeof(double));
    solpcg->SolverCg(p, b, 0);
    ass->Adjoint(x, p);

    ass->Forward(p, c);
    solpcg->SolverCg(v, p, 0);

    ass->Adjoint(s, v);
    for (int i = 0; i < n; i ++)
        s[i] = c[i] - s[i];

    /* check 
    this->GetRhsBC(b, c, x, v, s);
    ass->PrintVector(b, m, "init b");
    ass->PrintVector(c, n, "init c");*/

    double delx = 0.0, xsum = 0.0;
    double dels = 0.0, ssum = 0.0;
    for (int i = 0; i < n; i ++) {
	if (delx > x[i]) delx = x[i];
	if (dels > s[i]) dels = s[i];
    }
    delx = -1.5 * delx;
    dels = -1.5 * dels;

    for (int i = 0; i < n; i ++) {
	x[i] += delx;
	s[i] += dels;
	xsum += x[i];
	ssum += s[i];
    }

    double pxs = dot_product(x, s, n);
    delx = 0.5 * pxs / ssum;
    dels = 0.5 * pxs / xsum;

    for (int i = 0; i < n; i ++) {
	x[i] += delx;
	s[i] += dels;
    }
    delete [] p;
}


void SolverLin::ModifyFuncPara (double *x0, double *v0, double *s0) {
    ass->UpdateX(x0, v0, s0);
}

void SolverLin::UpdateDs (double *ds, double *dx) {
    for (int i = 0; i < n; i ++)
	ds[i] = -1.0 * ass->s0[i] * (1.0 + dx[i]);
}

void SolverLin::GetRhsPredictor (double *ru, double *rb, double *b, double *c) {
   memset(ru, 0, n * sizeof(double));
   memset(rb, 0, m * sizeof(double));
   ass->Adjoint(ru, ass->v0);
   for (int i = 0; i < n; i++) {
       ru[i] = ass->x0[i] * (c[i] - ru[i]);
   }
   ass->Forward(rb, ass->x0);
   for (int i = 0; i < m; i++) {
       rb[i] = b[i] - rb[i];
   }
}


void SolverLin::GetRhsCorrector (double *ru, double *rb, double *dr, double *b, double *c) {
   memset(ru, 0, n * sizeof(double));
   memset(rb, 0, m * sizeof(double));
   ass->Adjoint(ru, ass->v0);
   for (int i = 0; i < n; i++) {
       ru[i] = ass->x0[i] * (c[i] - ru[i]) + dr[i];
   }
   ass->Forward(rb, ass->x0);
   for (int i = 0; i < m; i++) {
       rb[i] = b[i] - rb[i];
   }
}


void SolverLin::GetRhsBC (double *b, double *c, double *x, double *v, double *s) {
   memset(b, 0, m * sizeof(double));
   memset(c, 0, n * sizeof(double));
   ass->Forward(b, x);
   ass->Adjoint(c, v);
   for (int i = 0; i < n; i ++)
       c[i] = c[i] + s[i];
}

void SolverLin::EstimateDr (double *dr, double *dx, double *ds) {
    double mu_aff = 0.0;
    double alpha_pri = 1.0;
    double alpha_dual= 1.0;
    double alpha;
    for (int i = 0; i < n; i ++) {
            if (dx[i] < 0) {
	       alpha = -ass->x0[i]/dx[i];
	       if (alpha_pri > alpha)
		  alpha_pri = alpha;
	    }
            if (ds[i] < 0) {
	       alpha = -ass->s0[i]/ds[i];
	       if (alpha_dual > alpha)
		  alpha_dual = alpha;
	    }
    }
    for (int i = 0; i < n; i ++) {
        mu_aff += (ass->x0[i] + alpha_pri * dx[i]) 
                * (ass->s0[i] + alpha_dual * ds[i]);
    }
    mu_aff = mu_aff / this->mu / n;
    double sigma = mu_aff * mu_aff * mu_aff;
    for (int i = 0; i < n; i++) {
        dr[i] = dx[i] * ds[i] - sigma * this->mu;
    }
    cout << "  sigma = " << sigma << ", mu = " << mu << endl;
}

void SolverLin::UpdateXVS (double *x1, double *v1, double *s1, double *dx, double *dv, double *ds) {
    double alpha_pri = 100.0;
    double alpha_dual= 100.0;
    double alpha;
    for (int i = 0; i < n; i ++) {
         if (dx[i] < 0) {
	    alpha = -ass->x0[i]/dx[i];
	    if (alpha_pri > alpha)
	       alpha_pri = alpha;
	 }
         if (ds[i] < 0) {
	    alpha = -ass->s0[i]/ds[i];
	    if (alpha_dual > alpha)
	       alpha_dual = alpha;
	 }
    }
    alpha_pri *= this->eta;
    alpha_dual *= this->eta;
    if (alpha_pri > 1.0) alpha_pri = 1.0;
    if (alpha_dual> 1.0) alpha_dual= 1.0;

    for (int i = 0; i < n; i ++) {
        x1[i] = ass->x0[i] + alpha_pri * dx[i];
        s1[i] = ass->s0[i] + alpha_dual * ds[i];
    }
    for (int i = 0; i < m; i ++) {
        v1[i] = ass->v0[i] + alpha_dual * dv[i];
    }
}

void SolverLin::InnerPoint (double *x, double *v, double *s, double *b, double *c) {

    double *dx = new double[n];
    double *dv = new double[m];
    double *ds = new double[n];

    double *dr = new double[n];
    double *ru = new double[n];
    double *rb = new double[m];

    int pcg_niter = 10;
    double pcg_eps = 1e-06;
    double pcg_wns = 0.01;

    int iter_out = 0;
    int iter_inn = 0;
    double r_perc = 1.0;
    double r_x0 = 1.0, r_x1;

    double mis0, mis1;
    double r1, r2;

    double minloss;
    double maxloss;

    // create an output file
    ofstream fp("output/kktlp_matrix_test.txt");

    // define a CG operator used in computing dx
    SolverPCG *solpcg = new SolverPCG (pcg_eps, pcg_wns, pcg_niter, ass);
    // this->TestPCG(solpcg); exit(0);

    // get initial value 
    this->GetInitPoint(x, v, s, b, c, solpcg);
    /* check 
    ass->PrintVector(x, n, "init x");
    ass->PrintVector(v, m, "init v");
    ass->PrintVector(s, n, "init s"); */
    for (int i = 0; i < n; i ++)
        fp << x[i] << " ";
    for (int i = 0; i < m; i ++)
        fp << v[i] << " ";
    for (int i = 0; i < n; i ++)
        fp << s[i] << " ";
    fp << endl;

    minloss = dot_product(x, c, n);
    maxloss = dot_product(v, b, m);

    for (int iter = 0; iter < this->niter; iter ++) {

	cout << endl << " ### iter = " << iter << "/" << niter << ", min loss =" << minloss << ", max loss =" << maxloss << endl;

	// update parameters
        this->ModifyFuncPara (x, v, s);

	if (iprint == 1){};

	// predictor
        this->GetRhsPredictor(ru, rb, b, c);
	if (iprint == 1) {
           ass->PrintVector(ru, n, "ru");
           ass->PrintVector(rb, m, "rb");
	}

        solpcg->ProjCG (dx, dv, ru, rb);
        this->UpdateDs (ds, dx);
	if (iprint == 1) {
           ass->PrintVector(dx, n, "dx_aff");
           ass->PrintVector(dv, m, "dv_aff");
           ass->PrintVector(ds, n, "ds_aff");
	}

	/*
        ass->KKTfunc(ru, rb, dx, dv);
	if (iprint == 1) {
           ass->PrintVector(ru, n, "1 ru");
           ass->PrintVector(rb, m, "1 rb");
	} exit(0);*/


	// preconditionging
        ass->Precon(dx, dv, ds, 0);

	// estimate dr from affine 
        this->EstimateDr (dr, dx, ds);
	if (iprint == 1)
           ass->PrintVector(dr, n, "xs_aff");

	// corrector
        this->GetRhsCorrector(ru, rb, dr, b, c);
	if (iprint == 1) {
           ass->PrintVector(ru, n, "ru");
           ass->PrintVector(rb, m, "rb");
	}

        solpcg->ProjCG (dx, dv, ru, rb);
        this->UpdateDs (ds, dx);
	if (iprint == 1) {
           ass->PrintVector(dx, n, "dx");
           ass->PrintVector(dv, m, "dv");
           ass->PrintVector(ds, n, "ds");
	}


	// preconditionging
        ass->Precon(dx, dv, ds, 0);
	if (iprint == 1) {
           ass->PrintVector(dx, n, "dx");
           ass->PrintVector(dv, m, "dv");
           ass->PrintVector(ds, n, "ds");
	}

	// update variables
        this->UpdateXVS (x, v, s, dx, dv, ds);
	if (iprint == 1) {
           ass->PrintVector(x, n, "update x");
           ass->PrintVector(v, m, "update v");
           ass->PrintVector(s, n, "update s");
	}
        for (int i = 0; i < n; i ++)
            fp << x[i] << " ";
        for (int i = 0; i < m; i ++)
            fp << v[i] << " ";
        for (int i = 0; i < n; i ++)
            fp << s[i] << " ";
        fp << endl;

        minloss = dot_product(x, c, n);
        maxloss = dot_product(v, b, m);

	/*
	iter_inn = 0;
	dr = 1.0;
	step /= beta;

        cout << " --- --- outer = " << iter_out << "/" << outer << ", inner = 0/" << inner << ", step = " << step << ", mis=" << mis0 << endl;
	exit(0);

	while (dr > 0 && iter_inn < inner) {
          iter_inn += 1;

          VectorAdd (x1, x0, dx, 1.0, step, nm);
	  //ass.Constraint(x1);
	  if (iprint == 1) {
	     cout << endl << " --- inner update x ---- " << endl;
             exaio.printArray(x1, nm);
	  }

          ass.residual_aff(b1, x1);
          mis1 = ass.misfit(x1);
          r_x1 = norm2(b1, nm);
          if (iprint == 1) {
             cout << endl << " --- update residual ---- " << endl;
             exaio.printArray(b1, nm);
          }

	  dr = r_x1 - (1.0 - alpha * step) * r_x0;
          cout << " --- --- outer = " << iter_out << "/" << outer << ", inner = " << iter_inn << "/" << inner << ", step = " << step << " , mis=" << mis1 << ", r1/r0 = " << r_x1 << "/" << r_x0 << endl;
          step = beta * step;

        }

	if (iter_out % 10 == 1 && iprint == 1) {
	   cout << " --- outer update x ---- " << endl;
           exaio.printArray(x1, nm);
	}

	// update ass.x0
	ass.UpdateX (x1);
	ass.maxres(&r1, &r2, b1);
        r_perc = abs(r_x1 / r_x0 - 1.0);
        //r_perc = abs(mis1 / mis0 - 1.0);
	cout << " --- outer = " << iter_out << "/" << outer << " mis=" << mis1 << " r_perc = " << r_perc << " res()=" << r1 << " res(Ax-b)" << r2 << endl;

*/
    }

    fp.close();

    cout << endl;
    cout << " --- final x ---- " << endl;
    //exaio.printArray(x1, nm);
    
    delete [] dx;
    delete [] dv;
    delete [] ds;
    delete [] dr;
    delete [] ru;
    delete [] rb;
}


/* Test Projected CG operator */
void SolverLin::TestPCG(SolverPCG *sol) {

    // set a solution
    double *x = new double[n];
    double *v = new double[m];
    for (int i = 0; i < n; i ++)
        x[i] = rand()/(RAND_MAX + 1.0);
    for (int i = 0; i < m; i ++)
        v[i] = rand()/(RAND_MAX + 1.0);

    // get RHS of equation
    double *yc = new double[n];
    double *yb = new double[m];
    ass->KKTfunc(yc, yb, x, v);
    //ass->PrintVector(yc, n, "ccc");
    //ass->PrintVector(yb, m, "bbb");
	    

    // check cg oporator
    sol->CheckCgsOp (0);
    sol->CheckCgsOp (1);


    // solve the problem
    double *x0 = new double[n];
    double *v0 = new double[m];
    sol->ProjCG(x0, v0, yc, yb);

    // print result
    cout << "index   predicted   true " << endl;
    for (int i = 0; i < n; i ++) {
        cout << " i, x = " << i << ", " << x0[i] << ", " << x[i] << endl; 
    }
    for (int i = 0; i < m; i ++) {
        cout << " i, v = " << i << ", " << v0[i] << ", " << v[i] << endl; 
    }

    delete [] x;
    delete [] v;
    delete [] x0;
    delete [] v0;
    delete [] yc;
    delete [] yb;
}
