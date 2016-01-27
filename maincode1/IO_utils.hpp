#ifndef IOUTILS_H
#define IOUTILS_H

#include <vector>

//*************************************************************************

// idiosyncratic input/ output - lots of overloading

void outputBinary_2d(char *fname,std::vector<std::vector<double>> data,int nt,double dx,double dt,
                     int zmin);
void outputBinary_2d(char *fname,std::vector<std::vector<double>> data,int nt,double dx,double dt,
                     int zmin, int iend);
void outputBinary_2d(char *fname,std::vector<std::vector<double>> data,int nt,double dx,double dt,
                     int zmin, int iend, int jend);
void outputascii(char *fname,std::vector<double> a);
void outputascii(char *fname,std::vector<int> a);

//#########################################################################

// inputs

void inputascii(char *fname,int zmax,std::vector<double> &a,std::vector<double> &b,int &n);

std::vector<double> inputBinary(char *fname, int &fileLen, double &firstval);
/*  PURPOSE: read in binary data in double precision
 input:		fname: filename for binary data
 fileLen : length of file
 output:		double precision vector of values
 */

std::vector<double> inputBinaryOld(char *fname, int &fileLen, double &firstval);

//**************************************************************************

std::vector<std::vector<double>> inputBinary2D(char *fname, int &dim1, int &dim2, int NZ, int NX);
//  PURPOSE: reads in binary data from filename fname
//			and loads it to a NX x NZ array

void loadBCBFields(std::vector<std::vector<double>> &pxbcb, std::vector<std::vector<double>> &pxbcbn, std::vector<std::vector<double>> &pybcb,
                   std::vector<std::vector<double>> &pybcbn, std::vector<std::vector<double>> &vxbcb, std::vector<std::vector<double>> &vxbcbn,
                   std::vector<std::vector<double>> &vybcb, std::vector<std::vector<double>> &vybcbn,int iefbc,int jebc,int ntload);
//  PURPOSE: load back boundary fields needed for re-starting the FDTD code

void loadFields(std::vector<std::vector<double>> &r, std::vector<std::vector<double>> &rn, std::vector<std::vector<double>> &vx, std::vector<std::vector<double>> &vxn,
                std::vector<std::vector<double>> &vy, std::vector<std::vector<double>> &vyn,int NX,int NZ,double dx,int ntload);
//  PURPOSE: load fields needed for re-starting the FDTD code

#endif