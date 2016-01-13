
// functions having to do with Acoustic Propagation

const double clambda = 300.0;

//#########################################################################
// setup function for linear code

void acoustSetup(double *Vxacterm,double *Vyacterm, double *Pacterm,
		double *pcylterm, double *advecterm, double *shearterm,   
		double *rho, double *c1d, double *w1d, double dt, double dx);

/* PURPOSE: to set up the acoustic constants for the time-stepping loops
these values are updated within the routine
  Vxacterm : the term to update Vx 
  Vyacterm : the term to update Vy 
  Pacterm  : the term to update the Pressure 
  pcylterm : the cylindrical Pressure term 
  advecterm : the advection term to use for p, vr, and vz (no shear yet) 
the next are inputs, unchanged within
  rho : 1D density profile 
  c1d : 1D sound speed
  w1d : 1D wind speeds
  dt, dx   : discretization in time and space
*/

// use this for nonlinear code, should overload this later
void acoust_setupN(double *dtdxorz, double *dtdxrhocsq, double *pcylterm, 
				   double *rho, double *c1d, double dt, double dx);

//#########################################################################
// Use this routine to do the first time step for linear code
// (not central differenced in time)

void TstepOne(double **vx, double **vy, double **p, double **vxn, 
			double **vyn, double **pn, double *Vxacterm, double *Vyacterm, 
			double *Pacterm, double *pcylterm, double *rvec, int *jztopo,
			int *jxtopo, double dt, double dx, int widx);

// Use this routine to do the first time step for nonlinear code

void Tstep_OneN(double **vx,double **vy,double **r,double **vxn,
				double **vyn, double **rn, double **pn, double *rho,  
				double *cc, double dt, double dx, int widx);

//#########################################################################

// for linear propagation use either 
// acoustTerm (2nd order accuracy) or acoustTerm4 (4th order)
// could overload acoustTerm s.t. it did the ENTIRE matrix by default
void acoustTerm(double **pacoust, double **vxacoust, double **vyacoust, 
				double **pn, double **vxn, double **vyn, double *dtdxorho,
				double *dtdxorz, double *dtdxrhocsq, double *pcylterm,
				double *wterm, double *swterm, double *rvec, int iwindoleft, 
				int iwindoright, int jwindotop, bool ileft, bool iright);

void acoustTerm4(double **pacoust, double **vxacoust, double **vyacoust, 
				double **pn, double **vxn, double **vyn, double *dtdxorho,
				double *dtdxorz, double *dtdxrhocsq, double *pcylterm,
				double *rvec, int iwindoleft, int iwindoright, bool ileft);

// use this one for Nonlinear propagation
void NonlinAcoust(double **rnacoustnl, double **vxacoustnl, 
				  double **vyacoustnl, double **pn, double **rn, double **vxn, 
				  double **vyn, double *rho, double *rvec,  double d2tdx, 
				  double d4tdx, double d2t, double dx, int iwindoleft, 
				  int iwindoright, bool ileft);

//#########################################################################  


