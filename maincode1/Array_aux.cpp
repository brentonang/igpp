#include <cmath>
#include <cstdio>          
#include <ctime>           // needed for timestamp
#include <iostream>
#include <iomanip> 
#include <stdlib.h>
 
#include "Array_aux.hpp"

using namespace std;

/*###############################################################################
 
            THE FIRST SEVERAL FUNCTIONS ARE PRETTY GENERALLY USEFUL
 
###############################################################################*/
// void printvector(double *a,int n) {
void printvector(vector<double> a) { // is n necessary? will there be cases where not the full vector size will be printed out?
    for (int i = 0; i <= a.size(); i++) printf("%16.12f\n",a[i]);
}

//#############################################################################*/
// void printvector(int *a,int n) {
void printvector(vector<int> a) {
    for (int i = 0; i <= a.size(); i++) printf("%d\n",a[i]);
}

/*###############################################################################
/##### find the minimum value of a vector in the range i1 through i2 ##########*/
/**********************************DEPRECATED*************************************/
// Use std::min_element, needs #include <algorithm> e.g. 
// vector<double> vector1 {4, 6, 8, 9, 11, 2, 66, 1000000, -2};
// cout << "The smallest element of vector1 is " << *std::min_element(vector1.begin(), vector1.end()) << endl;

// double vector_min(double *vector, int &imin,int i1,int i2)
// {
//     int i; imin = i1;
//     double vecmin=vector[i1];
//     for (i=i1+1; i<=i2; i++) {
//         if (vector[i]<vecmin ) {
//             vecmin = vector[i];
//             imin =i;
//         }
//     }
//     return vecmin;
// }

/*###############################################################################
/##### find the maximum value of a vector in the range i1 through i2 ##########*/
/**********************************DEPRECATED*************************************/
// Use std::max_element, needs #include <algorithm> e.g. 
// vector<double> vector1 {4, 6, 8, 9, 11, 2, 66, 1000000, -2};
// cout << "The largest element of vector1 is " << *std::max_element(vector1.begin(), vector1.end()) << endl;

// double vector_max(double *vector, int &imin,int i1,int i2)
// {
//     int i; imin = i1;
//     double vecmax=vector[i1];
//     for (i=i1+1; i<=i2; i++) {
//         if (vector[i]>vecmax ) {
//             vecmax = vector[i];
//             imin =i;
//         }
//     }
//     return vecmax;
// }

/*############################################################################### 
/######### dynamically allocate vector & initialize to 0 - overloaded #########*/
/***********************************DEPRECATED***********************************/
// Initialize vector using c++11 standard
// vector<double> vector1(n) where n = size
// automatically initializes a vector1 of size n to 0s e.g.
// vector<double> vector1(15) -> double vector of size 15 initialized to 0s
// vector<double> vector1(15, 5) -> double vector of size 15 initialized to 5s

// double * zeros(int n)
// {
//     double *vector;
// 	vector = initialize(n,0.0);
//     return (vector);
// }

/*-------------------------------------------------------------------------------
/--------- dynamically allocate array & initialize to 0 - overloaded ----------*/
/***********************************DEPRECATED***********************************/
// Initialize 2D vector using c++11 standard
// vector<vector<double>> vector1(row, vector<double>(col, value)); e.g.
// vector<vector<double>> vector1(15, vector<double>(14); -> 2D double vector 15x14 of 0s
// vector<vector<double>> vector1(22, vector<double>(12, 3)); -> 2D double vector 22x12 of 3s

// double ** zeros(int nx, int my)
// {
// 	double **array;
//     array = initialize(nx,my,0.0);
// 	return (array);
// }

/*############################################################################### 
/##### dynamically allocate array & initialize to given value - overloaded ####*/
/***********************************DEPRECATED***********************************/
// See above function explanation

// double ** initialize(int nx, int my, double value)
// {
//     int i,j;
// // printf("zeros2d: Allocating %lf MB memory\n",((nx+1)*(my+1))/((double)1e6));
// // setup array
//     double **arr=new double*[nx+1];
//     for(i=0; i<=nx; ++i ) arr[i] = new double[my+1];
//     // initialize, add dimension length in 0 elements
//     for (i=1; i<=nx; i++) for (j=1; j<=my; j++) arr[i][j]=value;
//     for (i=0; i<=nx; i++) arr[i][0]=(double)my;
//     for (j=1; j<=my; j++) arr[0][j]=(double)nx;
//     return (arr);
// }
 
/*-------------------------------------------------------------------------------
/--- dynamically allocate a vector & initialize to given value - overloaded ---*/
/***********************************DEPRECATED***********************************/
// See above function explanation

// double * initialize(int n,double value)
// {
//     double *vector=new double[n+1];
//     // initialize, add dimension length in 0 elements
//     for (int i=1; i<=n; i++) vector[i] = value;
//     vector[0]=(double)n;
//     return (vector);
// }

/*###############################################################################
/############# much like initialize but outputs an integer vector #############*/
/***********************************DEPRECATED***********************************/
// See above function explanation e.g.
// vector<int> vector1(15, 3) -> vector of ints of size 15 of 3s

// int * iinitial(int n, int value)
// {
//     int *vector=new int[n+1];
//     // initialize, add dimension length in 0 element
//     for (int i=1; i<=n; i++) vector[i] = value;
//     vector[0] = n;
//     return (vector);
// }

/*###############################################################################
/##### dynamically allocate integer VECTOR & initialize to 0 - overloaded #####*/
/***********************************DEPRECATED***********************************/
// See above function explanation e.g.
// vector<int> vector1(15) -> vector of ints of size 15 os 0s

// int * izeros(int n)
// {
//     int *vector = new int[n+1];
//     for (int i=1; i<=n; i++) vector[i] = 0;
//     vector[0] = n;
//     return (vector);
// }

/*-------------------------------------------------------------------------------
--------------- dynamically allocate integer ARRAY & initialize to 0  ---------*/
/***********************************DEPRECATED***********************************/
// See above function explanation

// int ** izeros(int nx, int my)
// {
//     int i,j;
//     int **arr=new int*[nx+1];
//     for(i=0; i<=nx; ++i ) arr[i] = new int[my+1];
// // initialize, add dimension length in 0 elements
//     for (i=1; i<=nx; i++) for (j=1; j<=my; j++) arr[i][j]=0;
//     for (i=0; i<=nx; i++) arr[i][0]=(double)my;
//     for (j=1; j<=my; j++) arr[0][j]=(double)nx;
//     return (arr);
// }

/*###############################################################################
/######### swap vectors, dimensions must be listed in 0 rows & columns ########*/
/***********************************DEPRECATED***********************************/
// Use vector::swap e.g.
// vector<vector<double>> vector1(15, vector<double>(5));
// vector<vector<double>> vector2(2, vector<double>(2, 4));
// vector1.swap(vector2);

// void swap2d(double **a, double **b)
// {
//     int n= a[0][1], m = a[1][0], i1=1, j1=1;
//     swap2d(a,b,i1,n,j1,m);
//     return;
// }

/*-------------------------------------------------------------------------------
 ------------ swap 2d vectors, section of 1st dimension & all of 2nd  ---------*/
/***********************************DEPRECATED***********************************/
// Use std::swap_ranges, both vectors must have same size e.g.
// vector<vector<double>> vector1(15, vector<double>(15,3));
// vector<vector<double>> vector2(15, vector<double>(15));
//swap_ranges(vector1.begin()+3, vector1.end(), vector2.begin());

// void swap2d(double **a, double **b, int i1, int i2)
// {
//     int j1=1, m = a[1][0]; 
//     swap2d(a,b,i1,i2,j1,m);
//     return;
// }

/*-------------------------------------------------------------------------------
------------swap vectors, dimensions must be listed in 0 rows & columns -------*/
/***********************************DEPRECATED***********************************/
// See above function explanation

// void swap2d(double **a, double **b,int i1,int i2,int j1,int j2)
// {
//     int i,j;
//     double dum;
//     //   printf("swap2d: dimensions n = %d  m = %d\n",n,m);
//     for (i=i1; i<=i2; i++) {
//         for (j=j1; j<=j2; j++) {
//             dum=a[i][j];
//             a[i][j]=b[i][j];
//             b[i][j]=dum;
//         }
//     }
//     return;
// }

/*############################################################################### 
/######### round_to_sigdigs - rounds *DOWN* to N significant digits ###########*/

double round_to_sigdigs(double value, int N)
{
    double factor = pow(10.0, N - ceil(log10(fabs(value))));
    return floor(value * factor) / factor; 

/* do the next few lines to round properly
    value *= factor ;
    value = value >= 0.0 ? floor(value + 0.5) : ceil(value - 0.5);
    return value / factor; */    
}

/*############################################################################### 
/####################### timestamp - returns date and time ####################*/
void timestamp()

// see http://www.cplusplus.com/reference/ctime
// also see p1262,1263 of Stroustrup book 
{
        
    time_t rawt2=time(NULL);
    tm *ptime = localtime(&rawt2);
    //    printf ("Current local time and date: %s", asctime(ptime));
    char buffer [80];
    strftime (buffer,80,"Local time: %D, %T \n",ptime);
    puts (buffer);
    return;
}

/*###############################################################################
/###############  linear_interp, assumes xin & xout are sorted ################*/
// XXX should probably setup an error message in case xin or xout are not sorted

double *linear_interp(double xin[],double yin[], double xout[])

{
    double dum;
    int iout, iin;
    int nin = xin[0];
    int nout = xout[0];
    double *yout = zeros(nout);
    int i1=1;
    
    if (nin!=yin[0]) {
        printf("ERROR: problem in linear_interp \n");
        printf("Input X and Y vectors are not equal in length, exiting \n");
        exit(-1);
    }
    
    for (iout = 1; iout <= nout; iout++ )    {
        if (xout[iout] < xin[1]){
            yout[iout] = yin[1];
        }
        else if (xout[iout]>xin[nin])
        {
            yout[iout] = yin[nout];
        }
        else  {
            for ( iin = i1; iin < nin; iin++ ) {
                if ( xout[iout] >= xin[iin] && xout[iout] < xin[iin+1])    {
                    yout[iout] = yin[iin] + (xout[iout]-xin[iin]) *
                    (yin[iin+1] - yin[iin])/(xin[iin+1] - xin[iin]);  //slope
                    break;
                }
            }
        }
        i1=iin;      // program will work for an unsorted xout without this line
    }
    return yout;
}

/*###############################################################################
/###############  cubic_spline, assumes xin & xout are sorted ################*/
// NOTE: should include an error message in case xin or xout are not sorted
// see Numerical Recipes for similar code, but not the same

double *cubic_spline(double xin[],double yin[], double xout[], double y1d,
					 double ynd)

{
    int nin = xin[0];
    int nout = xout[0];
    double *yout = zeros(nout), *y2 = zeros(nout), *u = zeros(nout);
	int iin, iout,k;
	double sig, p, qn, un;
	
	// check that input xin and yin are the same size
    if (nin!=yin[0]) {
        printf("Bad xin and yin input to routine linear_interp, exiting \n");
        exit(-1);
    }
	
    if (y1d > 0.99e30) {
		y2[1] = 0.0; u[1] = 0.0;
	} else {
		y2[1] = -0.5;
		u[1] = (3./(xin[2]-xin[1]))*(yin[2]-yin[1])/(xin[2]-xin[1]-y1d);
	}
	
	for (iin = 2; iin <= nin-1; iin++ )    {
		sig = (xin[iin]-xin[iin-1])/(xin[iin+1]-xin[iin-1]);
		p=sig*y2[iin-1]+2.0;
		y2[iin]=(sig-1.0)/p;
		u[iin] = (yin[iin+1]-yin[iin])/(xin[iin+1]-xin[iin]) 
			- (yin[iin]-yin[iin-1])/(xin[iin]-xin[iin-1]);
		u[iin] = (6.0*u[iin]/(xin[iin+1]-xin[iin-1])-sig*u[iin-1])/p;	
	}		
	
    if (ynd > 0.99e30) {
		un = 0.0; qn = 0.0;
	} else {
		qn = 0.5;
		un = (3./(xin[nin]-xin[nin-1]))*(ynd-(yin[nin]-yin[nin-1])/
										 (xin[nin]-xin[nin-1]));
	}		
	y2[nin]=(un-qn*u[nin-1])/(qn*y2[nin-1]+1.);			
	
    for (k = nin-1; k >= 1; k-- )    {
		y2[k]=y2[k]*y2[k+1]+u[k];
		//		printf("2nd derviative y2[%d] = %14.12f, xin[%d] = %lf\n",k,y2[k],k,xin[k]);
	}
	delete[] u;
	
	// now find yout values at each xout
	int klo,khi;
    double h,b,a; 
	
	klo = 1;			// put this line inside the loop if it's unsorted
    for (iout = 1; iout <= nout; iout++ )    {
        if (xout[iout] < xin[1]){		// if xout is out of bounds
            yout[iout] = yin[1];
        }
        else if (xout[iout]>xin[nin]) {
            yout[iout] = yin[nout];
        }
        else  {
			khi = nin;			
			while (khi-klo > 1) {		// find interval by bisection
				k=(khi+klo) >> 1;
				if (xin[k] > xout[iout]) khi=k;
				else klo=k;
			}
			
		    h=xin[khi]-xin[klo];
		    if (h==0.0) {
			    printf("Bad xin input to routine cubic_spline,exiting \n");
			    exit(-1);
		    }
			
		    a=(xin[khi]-xout[iout])/h;
		    b=(xout[iout]-xin[klo])/h;
		    yout[iout] = a*yin[klo]+b*yin[khi]+( (a*a*a-a)*y2[klo]+
												(b*b*b-b)*y2[khi] ) * (h*h) / 6.0;
		}
	}
    return yout; 
}

//###############################################################################
//###############################  reshape   ##################################*/
// reads in a 1D vector and transforms it into a 2D array
// http://www.ics.uci.edu/~dan/class/165/notes/memory.html) (memory allocation)

double ** reshape(double *c,int nx,int mz)
{
    int ii,jj,kk=0;
    double **a = new double*[mz+1];
    for (ii=0; ii<mz+1; ii++) {a[ii] = new double[nx+1];}  //allocate array space
    
    //    printf("reshape: %d %d matrix\n",nx,mz);
    //see http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
    // for how to delete this too (and has simpler ways of initializing)
    
    for (ii=0; ii<= mz; ii++) {
        for (jj=0; jj<=nx; jj++) {
            a[ii][jj] = c[kk];
            kk=kk+1;
        }
    }
    
    return (a);
}

//###############################################################################
//###############################  reshapeFill   ##############################*/
// reads in a 1D vector and transforms it into a 2D array
// it also fills the matrix will zeros if the input 1D vector is too small

double ** reshapeFill(double *c, int nx, int mz, int bigMz)
{
    int ii,jj,kk=0;
    double **a = new double*[bigMz+1];				//allocate array space
    for (ii=0; ii<bigMz+1; ii++) {a[ii] = new double[nx+1];}
    
    printf("reshapeFill: %d x %d matrix\n",bigMz+1,nx+1);
    printf("\t %d x %d filled\n",mz+1,nx+1);
    
    for (ii=0; ii<= mz; ii++) {
        for (jj=0; jj<=nx; jj++) {
            a[ii][jj] = c[kk];
            kk=kk+1;
        }
    }
    
    //	printf("filling \n");
    
    double dimy = (double)bigMz;
    for (ii=mz+1; ii<= bigMz; ii++) {
        a[ii][0] = dimy;
        for (jj=1; jj<=nx; jj++) {
            a[ii][jj] = 0.0;
        }
    }
    
    return (a);		
}
//############################################################################### 
