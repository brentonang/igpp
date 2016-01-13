#include<iomanip> 
#include <cmath>
#include <ctime>           // needed for timestamp
#include <cstdio>          // needed for outputBinary_2d
#include <fstream>
using namespace std;

# include "Array_aux.hpp"
# include "IO_utils.hpp"

/*######################################################################### 
 
 FUNCTIONS SPECIFIC TO THE WAVE PROPAGATION PROBLEM
 
###############################d#########################################*/

/*#########################################################################
/##########  write Cstyle binary output - double precision   ############*/

// the default is to output the entire file

void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt,int zmin)
{
    int nx = (int) data[0][1];
	int mz = (int) data[0][0];
	outputBinary_2d(fname,data,nt,dx,dt,zmin,nx,mz);
    return;
}

//-------------- overloaded, only output a section of the file ------------
//-------------------------------------------------------------------------


void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt,
					 int zmin, int iend)
{
    int mz = (int) data[0][0];
	outputBinary_2d(fname,data,nt,dx,dt,zmin,iend,mz);
    return;
}

//-------------- overloaded, only output a section of the file ------------
//-------------------------------------------------------------------------

// the default is to output the entire file
// see http://www.cplusplus.com/reference/cstdio/

void outputBinary_2d(char *fname,double **data,int nt,double dx,double dt, 
					 int zmin, int iend, int jend)
{
    int i, mz;
	FILE *fid;
    
    //  error handling
    fid=fopen(fname,"wb");
	if (!fid) {
        printf(" \n\n  in outputBinary_2d - FATAL ERROR!\n");
		printf("  Unable to open the output file\n\n");
        exit(-1);
	}
    
    // the 0 elements of the data are un-used in this program so the
    // header information can be embedded within the data[0][] row.
    //    printf(" dimensions iend = %d  jend = %d\n",iend,jend);
    data[0][2] = dx;
    data[0][3] = double(round(nt));
    data[0][4] = dt;
	data[0][5] = double(round(zmin));
	data[0][6] = double(round(iend));
	data[0][7] = double(round(jend));
    
    // write out a potion of the matrix and finish up
    for (i=0; i<=iend; i++) fwrite(&data[i][0],sizeof(double),jend+1,fid);
    fclose(fid);
    
    return;
}

/*#########################################################################
 /###############  outputascii - used for 1D profiles ####################*/

void outputascii(char *fname,double *a)
{
    FILE *fid;
    int n = (int)a[0];
    fid=fopen(fname,"w");
    if (fid!=NULL) {
		// write out the vector
        for (int i=1; i<=n; i++) fprintf(fid," %15.12f \n",a[i]);
        fclose(fid);
    }
    else {
        printf("in outputascii: Unable to open the output file\n\n");
	}
	
    return;
}

//!!!!!!!!!!!!!!!!!!  overloaded function  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void outputascii(char *fname,int *a)
{
    FILE *fid;
    int n = (int)a[0];
    fid=fopen(fname,"w");
    if (fid!=NULL) {
        // write out the vector
        for (int i=1; i<=n; i++) fprintf(fid," %d \n",a[i]);
        fclose(fid);
    }
    else {
        printf("in outputascii: Unable to open the output file\n\n");
    }
    return;
}


/*#########################################################################
 /###############  inputascii - used for 1D profiles #####################*/

// this inputs ascii files in a very C-like manner, not C++
// check http://www.cplusplus.com/forum/general/30918/ for improvements
// and look at filestreams.pdf

void inputascii(char *fname,int zmax,double *&a,double *&b,int &n)
{
    FILE *fid;
    fid=fopen(fname,"r");
    if (!fid) {
        printf(" \n\n  in inputascii - FATAL ERROR!\n");
        printf("Unable to open the input file\n\n");
        exit(-1);
    }
    
    int i=0;
    double tmp,tmp2;
    // first determine length of file to read in (read to just above zmax)
    //    if (fid!=NULL) while (!feof(fid)) {
    while (!feof(fid)) {
        fscanf(fid,"%lf %lf",&tmp,&tmp2);
        i++;
        //                printf("%d   %lf  %lf\n",i,tmp,tmp2);
        if (tmp>(double)zmax) break;
    }
    rewind(fid);
    if (tmp<zmax)
        n=i-1;
    else
        n=i;
    //    printf(" n = %d \n",n);
    
    
    // initialize and read in the vectors
    a = zeros(n);  b = zeros(n);
    for (i=1; i<=n; i++) {
        fscanf(fid,"%lf %lf",&a[i],&b[i]); // start at i=1
    }
    fclose(fid);
    
    return;
}

/*#########################################################################
 /#######################  loadBCBFields  ################################*/
void loadBCBFields(double **&pxbcb, double **&pxbcbn, double **&pybcb,
                   double **&pybcbn, double **&vxbcb, double **&vxbcbn,
                   double **&vybcb, double **&vybcbn,int iefbc,int jebc,int ntload)
{
    char infilename[80];
    int nxx, nyy;
    
    // read in fields for the back PML boundary conditions
    
    // read in px files
    sprintf(infilename,"pxbcbn%2.2d.bin",ntload);
    printf("\n reading field from file %s \n",infilename);
    pxbcbn = inputBinary2D(infilename, nyy, nxx, jebc,iefbc);
    
    sprintf(infilename,"pxbcb%2.2d.bin",ntload);
    printf("nt = %d, reading from file %s \n",ntload,infilename);
    pxbcb = inputBinary2D(infilename, nyy, nxx, jebc,iefbc);
    
    // read in py files
    sprintf(infilename,"pybcb%2.2d.bin",ntload);
    printf("nt = %d, reading from file %s \n",ntload,infilename);
    pybcb = inputBinary2D(infilename, nyy, nxx, jebc,iefbc);
    
    sprintf(infilename,"pybcbn%2.2d.bin",ntload);
    printf("\n reading field from file %s \n",infilename);
    pybcbn = inputBinary2D(infilename, nyy, nxx, jebc,iefbc);
    
    
    // read in radial velocity files
    sprintf(infilename,"vxbcb%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vxbcb = inputBinary2D(infilename, nyy, nxx, jebc,iefbc+1);
    
    sprintf(infilename,"vxbcbn%2.2d.bin",ntload);
    printf("\n reading field from file %s \n",infilename);
    vxbcbn = inputBinary2D(infilename, nyy, nxx, jebc,iefbc+1);
    
    // read in vertical velocity files
    sprintf(infilename,"vybcb%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vybcb = inputBinary2D(infilename, nyy, nxx, jebc+1,iefbc+1);
    
    sprintf(infilename,"vybcbn%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vybcbn = inputBinary2D(infilename, nyy, nxx, jebc+1,iefbc+1);
    
    return;
}

/*#########################################################################
 /#######################  loadFields  ###################################*/
void loadFields(double **&r, double **&rn, double **&vx, double **&vxn,
                double **&vy, double **&vyn,int NX,int NZ,double dx,int ntload)
{
    char infilename[80];
    int nxx, nyy;
    
    // read in resistivity files
    sprintf(infilename,"r%2.2d.bin",ntload);
    printf("nt = %d, reading from file %s \n",ntload,infilename);
    r = inputBinary2D(infilename, nyy, nxx, NZ,NX);
    //  printf("loadFields: NX = %d, NZ = %d, dx=%lf\n\n",nxx,nyy,r[0][2]);
    
    if (dx!=r[0][2]){
        printf("FATAL ERROR in loadFields: dx is incorrect, exit\n");
        printf("    file : %s\n",infilename);
        exit(-1);
    }
    
    sprintf(infilename,"rn%2.2d.bin",ntload);
    printf("\n reading field from file %s \n",infilename);
    rn = inputBinary2D(infilename, nyy, nxx, NZ,NX);
    
    // read in radial velocity files
    sprintf(infilename,"vx%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vx = inputBinary2D(infilename, nyy, nxx, NZ,NX+1);
    
    sprintf(infilename,"vxn%2.2d.bin",ntload);
    printf("\n reading field from file %s \n",infilename);
    vxn = inputBinary2D(infilename, nyy, nxx, NZ,NX+1);
    
    
    // read in vertical velocity files
    sprintf(infilename,"vy%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vy = inputBinary2D(infilename, nyy, nxx, NZ+1,NX);
    
    sprintf(infilename,"vyn%2.2d.bin",ntload);
    printf("\n reading from file %s \n",infilename);
    vyn = inputBinary2D(infilename, nyy, nxx, NZ+1,NX);
    
    return;
}

/*#########################################################################
 /#######################  inputBinary2D  ################################*/
double ** inputBinary2D(char *fname, int &dim1, int &dim2, int NZ, int NX)
{
    double *dumv;
    int fileLen;
    double firstval;
    double **array;
    
//    dumv=inputBinary(fname, fileLen,firstval);
    dumv=inputBinaryOld(fname, fileLen,firstval);
    dim1 = (int)firstval;
    dim2 = fileLen/(dim1+1) - 1;
    
    if ( dim1 != NZ | dim2>NX){
        printf("FATAL ERROR in inputBinary2D: dimension mismatch \n");
        printf("    file : %s\n",fname);
        printf("dim1 = %d, NZ = %d\n",dim1,NZ);
        printf("dim2 = %d, NX = %d\n",dim2,NX);
        exit(-1);
    }
    
    if ( dim2==NX){
        printf(" calling reshape\n");
        array=reshape(dumv,dim1,dim2);
    } else {
        printf(" calling reshapeFill\n");
        array=reshapeFill(dumv,dim1,dim2,NX);
    }
    
    /*
     int imax = fmin(7,dim2);
     printf("in inputBinary2D, imax = %d\n",imax);
     for (int i=0; i<=imax; i++) {
     for (int j=0; j<=5; j++) {
     printf("%lf  ",array[i][j]);
     }
     printf("\n");
     }
     exit(1);
     */
    
    //	delete []dumv;
    
    return (array);
}

/*#########################################################################
 /#######################  inputBinary   #################################*/
// this code is derived from http://www.cplusplus.com/doc/tutorial/files/
// & from
//stackoverflow.com/questions/23066176/read-double-type-data-from-binary-file

double * inputBinary(char *fname, int &fileLen, double &firstval)
{
    char *memblock;
    double *double_values;
    streampos size;
    ifstream file (fname, ios::in|ios::binary|ios::ate);
    if (file.is_open())
    {
        size = file.tellg();
        memblock = new char [size];
        file.seekg (0, ios::beg);
        file.read (memblock, size);
        file.close();
        
        double_values = (double*)memblock;  //reinterpret as doubles
        fileLen = (int)size/sizeof(double);
        //        cout << "file length in inputBinary: " << fileLen << endl;
        firstval = double_values[0];
        delete[] memblock;
    }
    else {
        printf("FATAL ERROR: Unable to open file in inputBinary\n\n");
        exit(-1);
    }
    
    return(double_values);
}

/*#########################################################################
 /#######################  inputBinaryOld   ##############################*/

double * inputBinaryOld(char *fname, int &fileLen, double &firstval)
{
    struct rec
    {
        double x;
    };
    
    int counter, recsize;
    FILE *fid;
    struct rec my_record;
    
    //............. read in all values ..................
    fid=fopen(fname,"rb");
    if (!fid) {
        printf("Unable to open file\n\n");
        //        return;
    }
    recsize = sizeof(struct rec);  // machine-dependent #bytes in a float
    
    fseek(fid, 0, SEEK_END);
    fileLen=(int)ftell(fid)/recsize;   // #  values
    rewind(fid);                //	fseek(fid, 0, SEEK_SET);  clearerr(fid);
    
    double *tmp = new double[fileLen];        //dynamically allocate tmp, delete it in calling routine
    
    for ( counter=0; counter <= fileLen-1; counter++) {
        fread(&my_record,sizeof(struct rec),1,fid);
        tmp[counter] = my_record.x;
    }   
    fclose(fid);
    firstval = tmp[0];
    
    return (tmp);
}
