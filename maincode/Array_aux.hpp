#ifndef ARRAYAUX_H
#define ARRAYAUX_H

//####################### function prototypes   #################################

void printvector(vector<double> a,int n);          //overloaded
void printvector(vector<int> a,int n);

/*  PURPOSE: to print out the first n vales of vector a

    USAGE: printvector(vector,20);
*/

//*******************************************************************************
double vector_min(vector<double> vector, int &imin,int i1,int i2);

/* PURPOSE: find the minimum value of a vector, between elements i1 through i2
    also returns the index at which the minimum value is found.
 
    USAGE: double ccmin = vector_min(ccv,imin,1,N);
 */

//*******************************************************************************
double vector_max(vector<double> vector, int &imin,int i1,int i2);

/* PURPOSE: find the maximum value of a vector, between elements i1 through i2
 also returns the index at which the maximum value is found.
 
 USAGE: double ccmax = vector_max(ccv,imin,1,N);
 */

//*******************************************************************************
vector<vector<double>> zeros(int n,int m);                 // overloaded function

/* PURPOSE: dynamically allocates an n+1 by m+1 matrix with matrix dimensions 
    placed on the 0th row and column
    the value m is placed on the matrix[i][0] column
    the value n is placed on the matrix[0][j] row
    All other elements are initialized to zero
 
    USAGE: double **array=zeros(NR,NZ);
*/
//-------------------------------------------------------------------------------
vector<double> zeros(int n);                          // overloaded function

/* PURPOSE: dynamically allocates an n+1 vector with the dimension n at the 
    vector[0] element. 
    All other elements are initialized to zero
 
    USAGE: double *vector=zeros(NR);
*/

//*******************************************************************************
vector<int> izeros(int n);                          // overloaded function

/* same as zeros(n) but for integers
 
 USAGE: int *vector=izeros(NR);
 */

//-------------------------------------------------------------------------------
vector<vector<int>> izeros(int nx, int my);
/* same as zeros(n,m) but for integers

 USAGE: int **array=izeros(NR,NZ);
 */

//*******************************************************************************
vector<vector<double>> initialize(int nx, int my,double value);      // overloaded function

/* PURPOSE: dynamically allocates an n+1 by m+1 matrix with matrix dimensions 
    placed on the 0th row and column
    the value m is placed on the matrix[i][0] column
    the value n is placed on the matrix[0][j] row
    All other elements are initialized to *value*
 
    USAGE: double *vector = initialize(NR,value);
*/

//-------------------------------------------------------------------------------
vector<double> initialize(int n, double value);               // overloaded function

/* PURPOSE: dynamically allocates an n+1 vector with the dimension n at the 
    vector[0] element. 
    All other elements are initialized to *value*
 
    USAGE: double **array = initialize(NR,NZ,value); 
*/

//*******************************************************************************

vector<int> iinitial(int n, int value);               // overloaded function

/* PURPOSE: dynamically allocates an integer n+1 vector with dimension n at the
 vector[0] element.
 All other elements are initialized to *value*
 
 USAGE: int **array = iinitial(NR,value);
 */
//*******************************************************************************
void swap2d(double **a, double **b);
void swap2d(double **a, double **b, int i1, int i2);                // overloaded
void swap2d(double **a, double **b,int i1,int i2,int j1,int j2);  // overloaded

/* PURPOSE: swaps 2 matrices
either
	USAGE:  
		swap2d(amatrix,bmatrix);
			swap the entire contents of 2d matrices amatrix and bmatrix
			(dimensions of matrix B must be >= those of A, not checked)
 or
		swap2d(amatrix,bmatrix,i1,i2,j1,j2;
			swap the contents of amatrix(i1:i2,:) & bmatrix(i1:i2,:)
 or
		swap2d(amatrix,bmatrix,i1,i2,j1,j2;
			swap the contents of amatrix(i1:i2,j1:j2) & bmatrix(i1:i2,j1:j2)
*/

//*******************************************************************************

double round_to_sigdigs(double value, int N);

/* PURPOSE: rounds the input value to N significant figures. ROUNDS DOWN!!!!!!!
    inputs:
        value: number to be rounded down to N significant digits
        N: number of significant digits wanted
    output:
        double precision value, rounded down to N significant digits
*/
//*******************************************************************************
void timestamp();

/* PURPOSE: outputs the local time in the format Local time: 11/06/13, 15:21:54
 
 USAGE: timestamp();
*/

//*******************************************************************************

double *linear_interp(double xin[],double yin[], double xout[]);
/*  PURPOSE: do linear interpolation given the input vectors xin and yin
 The function is linearly interpolated at locations given by the
 input vector xout.
 vectors xin, yin, xout all give dimensions in their [0]th element
 
 Note that both xin and xout MUST be sorted from smallest to largest.
 Also, if some xout values are out of bounds, the output values are
 set to either xin(1), if xout<xin[1], or xin[nin], if xout>xin[nout]

 */
//*******************************************************************************
double *cubic_spline(double xin[],double yin[], double xout[], double y1d,
					 double ynd);

/*  PURPOSE: do cubic spline interpolation given input vectors xin and yin
 The function is interpolated at locations given by the input vector xout.
 vectors xin, yin, xout all give dimensions in their [0]th element
 
 y1d, ynd : value of first derivative at first and last endpoints of the 
 interpolating function. If y1d or y1n >= 1e30, the routine sets the 
 corresponding boundary condition for a natural spline, with 0 2nd 
 derivatives at that boundary
 
 Note that both xin and xout MUST be sorted from smallest to largest.
 Also, if some xout values are out of bounds, the output values are
 set to either xin(1), if xout<xin[1], or xin[nin], if xout>xin[nout]
 
 */

//*******************************************************************************
double ** reshape(double *c,int n,int m);

/* PURPOSE: to reshape a vector into an array with dimensions n x m	. The input
 must have n*m elements. This is not checked within the function. */

double ** reshapeFill(double *c,int n,int m, int bigM);

#endif