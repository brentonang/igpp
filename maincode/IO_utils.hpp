//*************************************************************************

// idiosyncratic input/ output - lots of overloading

void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt,
                     int zmin);
void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt,
                     int zmin, int iend);
void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt,
                     int zmin, int iend, int jend);
void outputascii(char *fname,double *a);
void outputascii(char *fname,int *a);

//#########################################################################

// inputs

void inputascii(char *fname,int zmax,double *&a,double *&b,int &n);

double * inputBinary(char *fname, int &fileLen, double &firstval);
/*  PURPOSE: read in binary data in double precision
 input:		fname: filename for binary data
 fileLen : length of file
 output:		double precision vector of values
 */

double * inputBinaryOld(char *fname, int &fileLen, double &firstval);

//**************************************************************************

double ** inputBinary2D(char *fname, int &dim1, int &dim2, int NZ, int NX);
//  PURPOSE: reads in binary data from filename fname
//			and loads it to a NX x NZ array

void loadBCBFields(double **&pxbcb, double **&pxbcbn, double **&pybcb,
                   double **&pybcbn, double **&vxbcb, double **&vxbcbn,
                   double **&vybcb, double **&vybcbn,int iefbc,int jebc,int ntload);
//  PURPOSE: load back boundary fields needed for re-starting the FDTD code

void loadFields(double **&r, double **&rn, double **&vx, double **&vxn,
                double **&vy, double **&vyn,int NX,int NZ,double dx,int ntload);
//  PURPOSE: load fields needed for re-starting the FDTD code


