//#######################    global constants   ###########################

const double PI = 3.141592653589793, gammaa=1.40200000000, Rgas=287.05;

//#########################################################################
// these are specific to fdtd wave propagation

int init_source(int nr,int zmin,int nz,int srcz,double dx,double **a,
                double *c1d, float fsrc, int &js);

int newinit_source(int nr,int zmin,int nz,int srcz,double dx,double **a,
                double *c1d, float fsrc, int *jztopo, int &js);
//  PURPOSE: initialize the source given the spatial discretization

//-------------------------------------------------------------------------

void acoustic_halfspace(int rmax,int zmin,int zmax,float fsrc,int srcz,
                        float time_end);

void model_grids(double dx, int NR, int NZ, int zmin, double *rvec,
				 double *rvecx, double *zvec, double *zz);
//-------------------------------------------------------------------------
double compute_dx(double *ccv,double *wwv,int rmax,int zmin,int zmax,
                  float fsrc, int &NR,int &NZ, int ndx);
double compute_dx(double *ccv,double *wwv,int rmax,int zmin,int zmax,
                  float fsrc, int &NR,int &NZ);

// default ndx = 12 nodes per wavelength

//-------------------------------------------------------------------------
double compute_dt(double *cc,double *ww, double *rho,float Tlast,
				  float fsrc,double dx,double muvisc,int &NT);

double *compute_rho(double c1d[], int zmin,double dx, int NZ,double grav);
double *compute_rho(double c1d[], int zmin,double dx, int NZ, double grav,
					double pambient[]);

double *compute_rhoprof(double cin[], double zin[], int nin, double grav);
double *fixrhoprof(double zcc[], double rhoprof[], double rhomin);
void compute_pambient(double pamb[],double rho[],double cin[]);

//#########################################################################  


