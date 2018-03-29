/* flag = 1,  dz = P * dx         *
 * flag = 0,  dx = inv(P) * dz    */
/* Here P = inv(x0)               */
SolverLin::Precon(double *dx, double *dv, double *ds, int flag) {
    if (flag == 1) {
       for (int i = 0; i < n; i ++) {
	   dx[i] = dx[i] / ass->x0[i];
       }
    }
    else if (flag == 0) {
       for (int i = 0; i < n; i ++) {
	   dx[i] = dx[i] * ass->x0[i];
       }
    }
}
