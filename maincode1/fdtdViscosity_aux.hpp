//#########################################################################

// these functions are called only if viscosity > 0

//#########################################################################

void viscSetup(double *viscx,double *viscy,double *viscylr,double *viscylz, 
		double *viscylrr, double *muArt, double muvisc, double *rho, 
		double *c1d, double dt, double dx, int NR, int NZ, double fmax);

void viscSetupN(double *muArt, double *orvecxdx, double *orvecx2, 
        double *orvecdx, double *rho, double *c1d, double *rvec, 
        double *rvecx, double dt, double dx, double fmax, double muvisc);

void viscterm(double **vxvisc, double **vyvisc, double **vxn, double **vyn,
		double *viscx, double *viscy, double *viscylr, double *viscylz, 
		double *viscylrr, double *rvec, double *rvecx, int iwindoleft, 
		int iwindoright);

// could overload this so it would do the entire matrix by default
void viscxterm(double **vxvisc, double **vxn, double *viscx,  
		double *viscylr, double *viscylrr,  double *rvecx, int iwindoleft, 
		int iwindoright, int jwindotop);

void viscxtermN(double **vxviscN, double **vxn, double **rn, double *rho, 
        double *orvecxdx, double *orvecx2, double d2tmu, double odx2,
        int iwindoleft, int iwindoright);

// could overload this so it would do the entire matrix by default
void visczterm(double **vyvisc, double **vyn, double *viscy,  
		double *viscylz, double *rvec,  int iwindoleft, 
		int iwindoright, int jwindotop);

void viscztermN(double **vyviscN, double **vxn, double **rn, double *rho, 
        double *orvecdx, double d2tmu, double odx2, int iwindoleft, 
        int iwindoright);

void artifiviscP(double **pn, double **partvisc, double *scale,  
		int iwindoleft, int iwindoright);

void artifiviscVx(double **vxn, double **Vxartvisc, double *scale,  
		int iwindoleft, int iwindoright);

void artifiviscVy(double **vyn, double **Vyartvisc, double *scale,  
		int iwindoleft, int iwindoright);

void viscterm4(double **vxvisc, double **vyvisc, double **vxn,double **vyn,			  
		double *viscx, double *viscy, double *viscylr, double *viscylz, 
		double *viscylrr, double *rvec, double *rvecx, int iwindoleft, 
		int iwindoright);


//#########################################################################  


