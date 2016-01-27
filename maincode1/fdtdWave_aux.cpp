#include <cmath>
#include <cstdio>
#include <ctime>           // needed for timestamp
#include <iomanip> 
#include <stdlib.h>        // exit, EXIT_SUCCESS, EXIT_FAILURE

#include "Array_aux.hpp"
#include "fdtdWave_aux.hpp"

/*######################################################################### 
 
 FUNCTIONS SPECIFIC TO THE WAVE PROPAGATION PROBLEM
 
//#########################################################################
/############ compute_dx - for 1D environmental model ###################*/

double compute_dx(double *ccv,double *wwv,int rmax,int zmin,int zmax,
                  float fsrc,int &NR,int &NZ, int ndx)
{
	// find minimum wavelength, using minimum sound and wind speeds
    int icmin, iwmin;
    double ccmin = vector_min(ccv,icmin,1,(int)ccv[0]);
    double wwmin = vector_min(wwv,iwmin,1,(int)wwv[0]);
	//    printf("ccmin= %lf at %d,wwmin= %lf @ %d \n",ccmin,icmin,wwmin,iwmin);
//    double lambda = (ccmin+fmin(0,wwmin))/fsrc;
	double lambda = (ccmin+wwmin)/fsrc;
    double dx = round_to_sigdigs(lambda/ndx,2);
    NR = ceil(rmax/dx);
    NZ = ceil((zmax-zmin)/dx);
    
    return dx;
}

//!!!!!!!!!!!!!!!!!!  overloaded function  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double compute_dx(double *ccv,double *wwv,int rmax,int zmin,int zmax,
                  float fsrc,int &NR,int &NZ)
{
	// find minimum wavelength, using minimum sound and wind speeds
    int ndx = 12;
    double dx = compute_dx(ccv,wwv,rmax,zmin,zmax,fsrc,NR,NZ,ndx);
    
    return dx;
}

/*#########################################################################
/#############  model_grids - discretize the environmental model ########*/

void model_grids(double dx, int NR,int NZ,int zmin,double *rvec,
				 double *rvecx, double *zvec, double *zz)

/*      compute location vectors for P in radial and vertical directions
        interpolate wind and sound speed profiles to vertical grid
 rvec are cell center locations
 revecx are cell boundary locations
 zvec are cell center locations - vertical
 zz are cell boundaries - vertical
 */
{
    int i;
    for (i=1; i<=NR; i++) {
		rvec[i] = ((double)i-0.5)*dx;
		rvecx[i] = (double)(i-1)*dx;
	}
    for (i=1; i<=NZ; i++) zvec[i] = zmin+((double)i-0.5)*dx;
    for (i=1; i<=NZ+1; i++) zz[i] = zmin+((double)i-1) * dx;
    
    return;
}

/*#########################################################################
/############ compute_dt - for 1D environmental model ###################*/

double compute_dt(double *cc,double *ww,double *rho,float Tlast,float fsrc,
				  double dx,double muvisc,int &NT)
{
// find maximum static sound speed, taking attenuation into account
	
    int icmax, iwmax, i, icq;
    double ccmax, wwmax, dt, q, ccqmax, omega =2*PI*fsrc;
	int iend = (int)cc[0];
	
	double *ccq = zeros(iend);	
	for (i=1; i<=iend; i++) {
		q = omega*muvisc/(rho[i]*cc[i]*cc[i]);
		ccq[i]	= cc[i]*(1 +0.375*q*q);
		//		printf("%d: q = %lf, cc = %lf ccq = %lf \n",i,q,cc[i],ccq[i]);
	}
	ccqmax = vector_max(ccq,icq,1,(int)ccq[0]);
	printf("ccqmax = %lf at %d th point \n",ccqmax,icq);
	
//    ccmax = vector_max(cc,icmax,1,(int)cc[0]);
    wwmax = vector_max(ww,iwmax,1,(int)ww[0]);
	printf("wwmax = %lf at %d th point \n",wwmax,iwmax);
    
 //   dt = 0.5*dx/(wwmax+sqrt(wwmax*wwmax+3*ccqmax*ccqmax));
	dt = 0.5*dx/sqrt(3)/(ccqmax+wwmax);
    dt = round_to_sigdigs(dt,2);
    NT = ceil(Tlast/dt);
    
// compute the actual spatial minimum discretization, 
// could comment them out without affecting the program
    double  ceffmin, dnn;    
    double *ceff=zeros(iend);
    for (i=1; i<=iend; i++) ceff[i]=ccq[i]+ww[i];
    
    ceffmin = vector_min(ceff,icmax,1,(int)cc[0]);  
    dnn = ceffmin/fsrc/dx;
    printf("minimum gridpoints per wavelength = %lf \n\n",dnn);
	
	delete [] ccq;
    
    return dt;
}

/*#########################################################################
/####################  compute_rho   ####################################*/

double *compute_rho(double c1d[], int zmin, double dx, int NZ,double grav, 
					double pambient[])
// compute rho at cell centers
{
	int irlo, i, iflag = 0;
	double dum=0.0, rhomin=0.00000003;
	double *temp = zeros(NZ), *gravh = initialize(NZ,grav);
    double *rho = zeros(NZ), *a = zeros(NZ);
	double pamb = pambient[1], EarthRadius = 6371.0;
	char ofilename[80];
	
    temp[1] = c1d[1]*c1d[1]/(gammaa*Rgas);     
 	a[1] = (zmin+dx/2)*grav/(Rgas*temp[1]);
    rho[1]=pamb*exp(-a[1])/(Rgas*temp[1]);
    for (i=2; i<=NZ; i++) {
        temp[i] = c1d[i]*c1d[i]/(gammaa*Rgas);
        gravh[i] = gravh[i]*pow((EarthRadius/(EarthRadius+i*dx/1000)),2.0);
        dum = dx*gravh[i]/(Rgas*temp[i]);
        a[i] = dum+a[i-1];
        rho[i]=pamb*exp(-a[i])/(Rgas*temp[i]);
		if (rho[i]<rhomin)	iflag = 1;		// damping flag for low rho
    }
	
	printf("\n last gravh in compute_rho = %lf \n",gravh[NZ]);
	
	if (iflag == 1) {// damp if rho is very very low
		double *b = zeros(NZ);
		printf("limit uppermost rho by altering gravity in compute_rho\n"); 
// could combine next 2 lines for a findNearest method
		for (i=1; i<=NZ; i++) b[i] = fabs(rho[i]-rhomin);
		dum = vector_min(b,irlo,1,NZ);
		printf("minimum rho was %12.9f at %f\n", rho[irlo],irlo*dx);		
				
		int i1 = fmax(irlo-40,1);
		for (i=i1; i<=NZ; i++) {
			gravh[i] = gravh[i]*(.4-.6*tanh((double)(i-irlo)/8.0));
			dum = dx*gravh[i]/(Rgas*temp[i]);
			a[i] = dum+a[i-1];
			rho[i]=pamb*exp(-a[i])/(Rgas*temp[i]);
		}
		dum = vector_min(rho,irlo,1,NZ);
		printf("minimum rho is now %12.9f at %f\n", rho[irlo],irlo*dx);
		delete [] b;
		printf("\n uppermost gravh in compute_rho = %12.9f \n",gravh[NZ]);
	}
	
	printf("\n uppermost rho in compute_rho = %12.9f \n",rho[NZ]);

	delete[] a; delete[] temp;
    return rho;
}
//!!!!!!!!!!!!!!!!!! overloaded function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double *compute_rho(double c1d[], int zmin,double dx, int NZ,double grav)
// compute rho at cell centers
{
	double *rho, *pambient = initialize(NZ,101325);
	rho = compute_rho(c1d,zmin,dx,NZ,grav, pambient);
	return rho;
}

/*#########################################################################
/##############  compute the atmospheric density profile ################*/

double *compute_rhoprof(double cin[], double zin[], int nin, double grav)

// compute rho at altitudes given in zin, not necessary uniformly spaced
{
	int i;
	double dum=0.0, temperature;
	double pstandard = 101325, EarthRadius = 6371000, gravh = grav;		
	double *rho = zeros(nin), *aexp = zeros(nin);
	
	temperature = cin[1]*cin[1]/(gammaa*Rgas); 
	gravh = grav*pow((EarthRadius/(EarthRadius+zin[1])),2.0);
 	aexp[1] = zin[1]*gravh/(Rgas*temperature);
    rho[1]=pstandard*exp(-aexp[1])/(Rgas*temperature);	
    for (i=2; i<=nin; i++) {		
        temperature = cin[i]*cin[i]/(gammaa*Rgas);
        gravh = grav*pow((EarthRadius/(EarthRadius+zin[i])),2.0);
        dum = (zin[i]-zin[i-1])*gravh/(Rgas*temperature);
        aexp[i] = dum+aexp[i-1];
        rho[i]=pstandard*exp(-aexp[i])/(Rgas*temperature);
    }
	delete[] aexp; 
	
	//	printf("\n uppermost rho in compute_rhoprof = %12.9f \n",rho[nin]);
    return rho;
}

/*#########################################################################
/################  "fix" the atmospheric density profile ################*/

double *fixrhoprof(double zcc[], double rhoprof[], double rhomin)

// fix the density profile so that values don't fall far below rhomin
{
	 int i, nsave, nsave2, nin = rhoprof[0];
	 double *rho = zeros(nin);
	 
	 for (i=1; i<=nin; i++) {
		 if (rhoprof[i]>3.0*rhomin) nsave=i;
		 if (rhoprof[i]<rhomin) break;
	 }
	 nsave2 = i;
	 printf("nsave=%d, nsave2=%d \n",nsave,nsave2);
	 
	 if ((nsave+3) >=nin) { //rhomin near top of the profile, leave as is
		 printf("leave the rho profile as is \n");
		 return rhoprof;
	 }
	 else {  // making alterations to the density profile
		 int iint = nsave+3;
		 double zplus = 12000;	// additional altitude, in metres
		 double *zzi = zeros(iint), *rhoi = zeros(iint);
		 for (i=1; i<=nsave; i++) {		// copy profile for rho< rhomin
			 zzi[i]=zcc[i];
			 rhoi[i]=rhoprof[i];
		 }
		 rhoi[nsave+1]=1.04*rhomin; zzi[nsave+1]=zcc[nsave2];
	 //		printf("copied rho up to zzi[%d]= %lf,  %14.11f\n",nsave,zzi[nsave],rhoi[nsave]);
		 if (zzi[nsave+1]+zplus > zcc[nin-1]) {			
			 double ratio = 1.+(zcc[nin-1]-zzi[nsave+1])/zplus;
			 zzi[iint-1] = zcc[nin-1]; rhoi[iint-1]=ratio*rhomin;
			 zzi[iint] = zcc[nin]; rhoi[iint]=ratio*rhomin;
			 printf("<10 km left in profile, rhotop = %14.11f\n",rhoi[iint]);			
		 } else {		// find point 10 km above rhomin altitude
			double *b = zeros(iint);
			double dum;
			int i10;
			for (i=1; i<=nin; i++) b[i] = fabs(zcc[i]-(zzi[nsave+1]+zplus));
			dum = vector_min(b,i10,1,nin);
			zzi[iint-1] = zcc[i10]; rhoi[iint-1]=2*rhomin;
			zzi[iint] = zcc[nin]; rhoi[iint]=2*rhomin;
			printf("zcc[%d] = %lf, zcc[%d] = %lf \n", nsave2, zcc[nsave2], i10, zcc[i10]);	
		 }
		 rho = cubic_spline(zzi,rhoi,zcc,0,0); 
		 return rho;
	 }
}

/*#########################################################################
/##############  compute the atmospheric pressure profile ###############*/

void compute_pambient(double pamb[],double rho[],double cin[])

// compute ambient pressure at nodes
{
	int nin = rho[0];
		
    for (int i=1; i<=nin; i++) {
        pamb[i]=rho[i]*cin[i]*cin[i]/gammaa;
    }
    return;
}

/*######################################################################### 
/########################## init_source - ###############################*/

int init_source(int nr,int zmin,int nz,int srcz,double dx,double **a,
                    double *c1d, float fsrc, int &js)
{
    int i,j;
	double rsource = 0.0;			// source range, m
    js = ceil(1+round((srcz-zmin)/dx - 0.5));
    double width=c1d[js]/fsrc;			// wavelength
    int na=4;
    int widx = ceil(na*width/dx);  
	printf("source wavlength = %lf, width = %d \n",width, widx);
    int is = floor(rsource/dx)+1;    
    int nr1 = widx + is;
    
    int jbegin = js-widx;		// int jbegin = js-widx< 1 ? 1: js-widx;
    int jend = js+widx;
	int ibegin = is-widx< 1 ? 1: is-widx;

    printf("ibegin = %d, jbegin = % d, jend = %d \n",ibegin,jbegin,jend);
    
// ERROR MESSAGES .........................................................
    if (jbegin<1) {
        printf(" \n\n in init_source - FATAL ERROR!\n");
        printf(" source altitude +width < minimum model altitude \n");
        printf(" srcalt-width = %lf, zmin = %d\n",zmin+(jbegin-1)*dx,zmin);
		printf(" HINT: increase srcalt or increase fsrc & recompile \n");
        exit(-1);
    }
    if ((nr1+1) > nr) {
        printf("\n\n in init_source - FATAL ERROR!\n");
        printf("source is at  %lf > max model range %lf\n",widx*dx,nr*dx);
		printf("HINT: increase rmax or increase frequency & recompile\n");
        exit(-1);
	}
    if (jend > nz) {
        printf(" \n\n in init_source - FATAL ERROR!\n");
        printf(" source altitude + width  > maximum model altitude \n");
        printf(" srcalt+width=%lf, zmax=%lf\n",zmin+(jend-1)*dx,zmin+(nz-1)*dx);
		printf(" HINT: increase zmax or decrease srcalt or increase fsrc \n");
		printf(" then recompile \n");
        exit(-1);
	}
// END ERROR MESSAGES .....................................................
    
// now add a (windowed) Gaussian source to the input vector. There is a 
// half-cell offset in i as source is on a cell boundary, p is @ center
    double R2 = width*width;
    double Rna = na*width, Rna2 = Rna*Rna;
    for (i=ibegin; i<=nr1; i++) {
        for (j=jbegin; j<=jend; j++) {
            double rad = dx*sqrt((j-js)*(j-js)+(i-is+0.5)*(i-is+0.5));
//            a[i][j] += rad > Rna ? 0: exp(-rad*rad/R2)*pow( (1-rad*rad/Rna2), 2);
            a[i][j] += rad > Rna ? 0: exp(-rad*rad/R2)*pow( cos(rad*PI/(2*Rna)), 2);
        }
    }
    
    return (nr1+1);      // compute 1st time step to this value, not to NR
}

/*######################################################################### 
/######################## newinit_source - ##############################*/

int newinit_source(int nr,int zmin,int nz,int srcz,double dx,double **a,
				double *c1d, float fsrc, int *jztopo, int &js)
{
    int i,j,jadd;
	double rsource = 0.0;			// source range, m
    js = fmax(ceil(1+round((srcz-zmin)/dx - 0.5)),1);
    double width=c1d[js]/fsrc;			// wavelength
    int na=4;
    int widx = ceil(na*width/dx);  
	printf("source wavlength = %lf, width = %d \n",width, widx);
    int is = floor(rsource/dx)+1;    
    int nr1 = widx + is;
    
    int jbegin = js-widx;		// int jbegin = js-widx< 1 ? 1: js-widx;
    int jend = js+widx;
	int ibegin = is-widx< 1 ? 1: is-widx;
    
// make sure that topography is flat in vicinity of the source
    int iflag = 0;
    for (i=ibegin+1; i<nr1; i++){
        if (jztopo[i] != jztopo[i-1]) iflag=1;
    }
	
    printf("ibegin = %d, jbegin = % d, jend = %d \n",ibegin,jbegin,jend);
	printf("source index in altitude, js = %d\n",js);
    
// ERROR MESSAGES .........................................................

    if (iflag != 0) {
        printf("\n\n in init_source - FATAL ERROR!\n");
        printf("topography must be flat from %lf to %lf m\n",(ibegin-1)*dx,nr1*dx);
        exit(-1);
    }
    
    if ((nr1+1) > nr) {
        printf("\n\n in init_source - FATAL ERROR!\n");
        printf("source is at  %lf > max model range %lf\n",widx*dx,nr*dx);
		printf("HINT: increase rmax or increase frequency & recompile\n");
        exit(-1);
	}
    if (jend > nz) {
        printf(" \n\n in init_source - FATAL ERROR!\n");
        printf(" source altitude + width  > maximum model altitude \n");
        printf(" srcalt+width=%lf \n",zmin+(jend-1)*dx);
		printf(" HINT: increase zmax or decrease srcalt & recompile\n");
        exit(-1);
	}
// END ERROR MESSAGES .....................................................
    
	// now add a (windowed) Gaussian source to the input vector. There is a 
	// half-cell offset in i as source is on a cell boundary, p is @ center
    double R2 = width*width;
    double Rna = na*width, Rna2 = Rna*Rna, rad, zdist,rdist;
	for (j=jbegin; j<=jend; j++) {
//		if (j>0) {          //assumes jztopo[1]=1
        if (j>jztopo[1]-1) {
			jadd = j;
		} else {
//			jadd = 1-j;     //assumes jztopo[1]=1
            jadd = 2*jztopo[1]-1-j;
		}
//		printf("jadd = %d\n",jadd);
		for (i=ibegin; i<=nr1; i++) {
			zdist = (j-js)*dx;
			rdist = i*dx-rsource;
			rad = sqrt(zdist*zdist+rdist*rdist);
//            a[i][jadd] += rad > Rna ? 0: exp(-rad*rad/R2)*pow( (1-rad*rad/Rna2), 2);
			a[i][jadd] += rad > Rna ? 0: exp(-rad*rad/R2)*pow( cos(rad*PI/(2*Rna)), 2);						
        }
    }
    
    return (nr1+1);      // compute 1st time step to this value, not to NR
}

//######################################################################### 
