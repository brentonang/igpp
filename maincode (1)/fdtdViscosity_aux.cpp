#include<iomanip> 
#include <cmath>
#include <ctime>           // needed for timestamp
#include <cstdio>          

# include "Array_aux.hpp"
# include "fdtdViscosity_aux.hpp"
# include "fdtdWave_aux.hpp"

/*######################################################################### 
 
 FUNCTIONS SPECIFIC TO WAVE PROPAGATION WITH VISCOSITY
 
###############################d#########################################*/
/*#########################################################################   
/## viscSetup - set constants related to viscosity for time step loop ###*/

void viscSetup(double *viscx,double *viscy,double *viscylr,double *viscylz, 
        double *viscylrr, double *muArt, double muvisc, double *rho, 
        double *c1d, double dt, double dx, int NR, int NZ, double fmax)

{	
// need to set up a viscylr and viscylz separately (1/rho)
    int i,j;
    double dum, dum2, omega = 2*PI*fmax, maxArtVisc=0.03;
	
	double c1 = 2*dt*muvisc/(dx*dx), c2 = dt*muvisc/dx;
	double c3 = 0.5*dt*muvisc/(omega*omega*dx*dx*dx*dx);

    // need to do checks on artVisc to ensure it's << 1
    for (j=1; j<=NZ; j++) {
        viscx[j] = c1/rho[j];					// atten term (2D)
     	viscylr[j] = c2/rho[j];					// visc cyl  (needs 1/r)
		muArt[j] = c3*c1d[j]*c1d[j]/rho[j];		// artificial atten
		if (muArt[j]>maxArtVisc) muArt[j] = maxArtVisc;
    }
	viscy[1] = viscx[1]; viscy[NZ+1] = viscx[NZ];
	viscylz[1] = viscylr[1]; viscylz[NZ+1] = viscylr[NZ];
	for (j=2; j<=NZ; j++) {
		viscy[j] = c1/((rho[j-1]+rho[j])/2);
		viscylz[j] = c2/((rho[j-1]+rho[j])/2);    // needs 1/r
	}
	
	viscylrr[1] = 0.0;  
	for (i=2; i<=NR; i++) viscylrr[i]= 2*dx/((double)(i-1)*dx);
	
//	printf("viscylr \n"); printvector(viscylr,15);printf("\n"); exit(1);
	
// getting maxima for viscylr and viscx
	dum = vector_max(viscx,i,1,NZ);
	printf("maximum viscx = %g at %d in viscSetup\n",dum,i);
	dum2 = vector_max(viscylr,j,1,NZ);
	printf("maximum viscylr = %g at %d in viscSetup\n",dum2,j);
	dum = vector_max(muArt,i,1,NZ);
	printf("maximum muArt = %g at %d in viscSetup\n\n",dum,i);
//	exit(1);
    return;
}

/*#########################################################################   
/## viscSetup - set constants related to viscosity for time step loop ###*/

void viscSetupN(double *muArt, double *orvecxdx, double *orvecx2, 
    double *orvecdx, double *rho, double *c1d, double *rvec, 
    double *rvecx, double dt, double dx, double fmax, double muvisc)

{	
    int i,j;
	int NZ = rho[0];
	int NR = rvec[0];
    double omega = 2*PI*fmax, maxArtVisc=0.03;	
	double c3 = 0.5*dt*muvisc/(omega*omega*dx*dx*dx*dx);
	
    // need to do checks on artVisc to ensure it's << 1
    for (j=1; j<=NZ; j++) {
		muArt[j] = c3*c1d[j]*c1d[j]/rho[j];		// artificial atten
		if (muArt[j]>maxArtVisc) muArt[j] = maxArtVisc;
    }
	
    orvecxdx[1] = 0; orvecx2[1] = 0;
	for (i=2; i<=NR; i++) {
        orvecxdx[i] = 1.0/(2*rvecx[i]*dx);
        orvecdx[i] = 1.0/(2*rvec[i]*dx);
        orvecx2[i] = 1.0/(rvecx[i]*rvecx[i]);
    }

    return;
}

/*######################################################################### 
/################################ viscterm ##############################*/

void viscterm(double **vxvisc, double **vyvisc, double **vxn, double **vyn,
		double *viscx, double *viscy, double *viscylr, double *viscylz, 
		double *viscylrr, double *rvec, double *rvecx, int iwindoleft, 
		int iwindoright)
{
    int i,j;
    int MR = (int) vyn[0][1];
    int MZ = (int) vxn[0][0];
//    printf("dimensions MR = %d MZ = %d \n",MR,MZ);
	
// Vx viscosity - vxvisc[1][j] = 0
// Vx is symmetric about z=0, anti-symmetric about r=0
    for (i=iwindoleft; i<=iwindoright; i++) {
        vxvisc[i][1] = viscx[1]*(vxn[i-1][1]+vxn[i+1][1]+vxn[i][2]
					- 3*vxn[i][1]) + viscylr[1]*(vxn[i+1][1]-vxn[i-1][1]
					- vxn[i][1]*viscylrr[i])/rvecx[i];
        for (j=2; j<=MZ-1; j++) { 
            vxvisc[i][j] = viscx[j]*(vxn[i-1][j]+vxn[i+1][j]
					+ vxn[i][j-1]+vxn[i][j+1]-4*vxn[i][j])
					+ viscylr[j]*(vxn[i+1][j]-vxn[i-1][j]
					- vxn[i][j]*viscylrr[i])/rvecx[i];
        }
        vxvisc[i][MZ] = viscx[MZ]*(vxn[i-1][MZ] + vxn[i+1][MZ]
					+ vxn[i][MZ-1] - 3*vxn[i][MZ])
					+ viscylr[MZ]*(vxn[i+1][MZ]-vxn[i-1][MZ]
					- vxn[i][MZ]*viscylrr[i])/rvecx[i];
    }
	
// Vy viscosity - vy is symmetric about r=0, antisymmetric about z=0
    for (i=iwindoleft; i<=iwindoright-1; i++) {
        for (j=2; j<=MZ; j++) { 
            vyvisc[i][j] = viscy[j]*(vyn[i-1][j]+vyn[i+1][j]
						+ vyn[i][j-1]+vyn[i][j+1]-4*vyn[i][j])
						+ viscylz[j]*(vyn[i+1][j]-vyn[i-1][j])/rvec[i];
        }
    }
	
    if (iwindoleft==2){           // if ileft==true
        for (j=2; j<=MZ; j++) {
            vyvisc[1][j] = viscy[j]*(vyn[1][j-1]+vyn[1][j+1]+vyn[2][j]
						- 3*vyn[1][j])      // vyn[1][j]= vyn[-1][j]
						+ viscylz[j]*(vyn[2][j]-vyn[1][j])/rvec[1];
        }		// skip *2 factor because of symmetry about i=1
    }	
    if (iwindoright==MR){           // if iright==true
        for (j=2; j<=MZ; j++) {
            vyvisc[MR][j] = viscy[j]*(vyn[MR][j-1] + vyn[MR][j+1]
						+ vyn[MR-1][j] - 3*vyn[MR][j]) + 2*viscylz[j]
						*(vyn[MR][j]-vyn[MR-1][j])/rvec[MR];
        }
    }
	return;
}

/*######################################################################### 
/################################ viscxterm #############################*/

void viscxterm(double **vxvisc, double **vxn, double *viscx, 
		double *viscylr, double *viscylrr,  double *rvecx, int iwindoleft, 
		int iwindoright, int jwindotop)
{
    int i,j;
    int MR = (int) vxn[0][1]-1;
    int MZ = (int) vxn[0][0];
	
//    printf("dimensions MR = %d MZ = %d \n",MR,MZ);
	
// Vx viscosity - vxvisc[1][j] = 0
// Vx is symmetric about z=0, anti-symmetric about r=0
    for (i=iwindoleft; i<=iwindoright; i++) {
        vxvisc[i][1] = viscx[1]*(vxn[i-1][1]+vxn[i+1][1]+vxn[i][2]
					- 3*vxn[i][1]) + viscylr[1]*(vxn[i+1][1]-vxn[i-1][1]
					- vxn[i][1]*viscylrr[i])/rvecx[i];
		if (jwindotop==MZ) {
			for (j=2; j<=MZ-1; j++) { 
				vxvisc[i][j] = viscx[j]*(vxn[i-1][j]+vxn[i+1][j]
					+ vxn[i][j-1]+vxn[i][j+1]-4*vxn[i][j])
					+ viscylr[j]*(vxn[i+1][j]-vxn[i-1][j]
					- vxn[i][j]*viscylrr[i])/rvecx[i];
			}
			vxvisc[i][MZ] = viscx[MZ]*(vxn[i-1][MZ] + vxn[i+1][MZ]
					+ vxn[i][MZ-1] - 3*vxn[i][MZ])
					+ viscylr[MZ]*(vxn[i+1][MZ]-vxn[i-1][MZ]
					- vxn[i][MZ]*viscylrr[i])/rvecx[i];
		} else {
			for (j=2; j<=jwindotop; j++) { 
				vxvisc[i][j] = viscx[j]*(vxn[i-1][j]+vxn[i+1][j]
					+ vxn[i][j-1]+vxn[i][j+1]-4*vxn[i][j])
					+ viscylr[j]*(vxn[i+1][j]-vxn[i-1][j]
					- vxn[i][j]*viscylrr[i])/rvecx[i];
			}
		}
	}
	return;
}

/*######################################################################### 
/############ viscxtermN  - nonlinear equivalent to viscxterm ###########*/

void viscxtermN(double **vxviscN, double **vxn, double **rn, double *rho, 
        double *orvecxdx, double *orvecx2, double d2tmu, double odx2,
        int iwindoleft, int iwindoright)

{
    int i,j;
    int MR = (int) vxn[0][1]-1;
    int MZ = (int) vxn[0][0];	
    double aa, bb;
//    printf("dimensions MR = %d MZ = %d in viscxtermN\n",MR,MZ);
	
// Vx is symmetric about z=0, anti-symmetric about r=0 - vxvisc[1][j] = 0		
    for (i=iwindoleft; i<=iwindoright; i++) {
        aa = orvecxdx[i];
        bb = orvecx2[i];
        
        vxviscN[i][1] = d2tmu/(rho[1]+rn[i][1]) * (odx2*
                (vxn[i-1][1] + vxn[i+1][1] +vxn[i][2]- 3*vxn[i][1])
                + aa*(vxn[i+1][1]-vxn[i-1][1]) - bb*vxn[i][1]);		
        for (j=2; j<=MZ-1; j++) { 
            vxviscN[i][j] = d2tmu/(rho[j]+rn[i][j]) * (odx2*(vxn[i-1][j] 
                + vxn[i+1][j] + vxn[i][j-1] +vxn[i][j+1]-4*vxn[i][j]) 
                + aa*(vxn[i+1][j]-vxn[i-1][j]) - bb*vxn[i][j]);
        }
		vxviscN[i][MZ] = d2tmu/(rho[MZ]+rn[i][MZ]) * (odx2*(vxn[i-1][MZ] 
                +  vxn[i+1][MZ] + vxn[i][MZ-1] - 3*vxn[i][MZ]) 
                + aa*(vxn[i+1][MZ]-vxn[i-1][MZ]) - bb*vxn[i][MZ]);
    }	
	return;
}

/*######################################################################### 
/############################### viscyterm ##############################*/

void visczterm(double **vyvisc, double **vyn, double *viscy,  
			   double *viscylz, double *rvec,  int iwindoleft, 
			  int iwindoright, int jwindotop)
{
    int i,j;
    int MR = (int) vyn[0][1];
    int MZ = (int) vyn[0][0]-1;
	int iend = int(fmin(MR-1,iwindoright));      
	
// Vy viscosity - vy is symmetric about r=0, antisymmetric about z=0
    for (i=iwindoleft; i<=iend; i++) {
        for (j=2; j<=jwindotop; j++) { 
            vyvisc[i][j] = viscy[j]*(vyn[i-1][j]+vyn[i+1][j]
						+ vyn[i][j-1]+vyn[i][j+1]-4*vyn[i][j])
						+ viscylz[j]*(vyn[i+1][j]-vyn[i-1][j])/rvec[i];
        }
    }	
    if (iwindoleft==2){           // if ileft==true
        for (j=2; j<=jwindotop; j++) {
            vyvisc[1][j] = viscy[j]*(vyn[1][j-1]+vyn[1][j+1]+vyn[2][j]
						- 3*vyn[1][j])      // vyn[1][j]= vyn[-1][j]
						+ viscylz[j]*(vyn[2][j]-vyn[1][j])/rvec[1];
        }		
    }	
    if (iwindoright==MR){           // if iright==true
        for (j=2; j<=jwindotop; j++) {
            vyvisc[MR][j] = viscy[j]*(vyn[MR][j-1] + vyn[MR][j+1]
						+ vyn[MR-1][j] - 3*vyn[MR][j]) + 2*viscylz[j]
						*(vyn[MR][j]-vyn[MR-1][j])/rvec[MR];
        }
    }
	return;
}

/*######################################################################### 
/############ viscztermN  - nonlinear equivalent to visczterm ###########*/


void viscztermN(double **vyviscN, double **vyn, double **rn, double *rho, 
            double *orvecdx, double d2tmu, double odx2, int iwindoleft, 
            int iwindoright)

{
    int i,j;
    int MR = (int) vyn[0][1];
    int MZ = (int) vyn[0][0]-1;	
    double aa;
//    printf("dimensions MR = %d MZ = %d in viscztermN\n",MR,MZ);
	int iend = int(fmin(MR-1,iwindoright));
	
// Vy viscosity - vy is symmetric about r=0, antisymmetric about z=0	
 
    for (i=iwindoleft; i<=iend; i++) {
        aa = orvecdx[i];
        for (j=2; j<=MZ; j++) { 
            vyviscN[i][j] = d2tmu/(rho[j]+rn[i][j]) * (odx2* (vyn[i-1][j] 
                    + vyn[i+1][j] + vyn[i][j-1]+vyn[i][j+1]-4*vyn[i][j])
                    + aa*(vyn[i+1][j]-vyn[i-1][j]));
        }
    }  
	
    if (iwindoleft==2){           // if ileft==true
        aa = orvecdx[1];
        for (j=2; j<=MZ; j++) {
            vyviscN[1][j] = d2tmu/(rho[j]+rn[1][j]) * (odx2*  
                    (vyn[2][j] + vyn[1][j-1]+vyn[1][j+1]-3*vyn[1][j])
                    + aa*(vyn[2][j]-vyn[1][j]));
        }
    }
    if (iwindoright==MR){           // if iright==true
       aa = orvecdx[MR];
       for (j=2; j<=MZ; j++) {
           vyviscN[MR][j] = d2tmu/(rho[j]+rn[MR][j]) * (odx2*  
                (vyn[MR-1][j] + vyn[MR][j-1]+vyn[MR][j+1]-3*vyn[MR][j])
                + aa*(vyn[MR][j]-vyn[MR-1][j]));
        }
    }
    
	return;
}

/*######################################################################### 
/############################## artifiviscP #############################*/

void artifiviscP(double **pn, double **partvisc, double *scale,  
			   int iwindoleft, int iwindoright)
// use symmetry of P about rigid surface and r=0
{
    int i,j;
    int MR = (int) pn[0][1];
    int MZ = (int) pn[0][0];
//    printf("dimensions MR = %d MZ = %d \n",MR,MZ);

	int ibegin = int(fmax(3,iwindoleft));
    int iend = int(fmin(MR-2,iwindoright));
	
// i = 1 terms
    if (iwindoleft<3) {
        partvisc[1][1] = scale[1]*(12*pn[1][1]
			-4*(pn[1][1]+pn[2][1]+pn[1][1]+pn[1][2])
			+(pn[2][1]+pn[3][1]+pn[1][2]+pn[1][3]));					
        partvisc[1][2] = scale[2]*(12*pn[1][2]
			-4*(pn[1][2]+pn[2][2]+pn[1][1]+pn[1][3])
			+(pn[2][2]+pn[3][2]+pn[1][1]+pn[1][4]));
        for (j=3; j<=MZ-2; j++) { 
            partvisc[1][j] = scale[j]*(12*pn[1][j]
                -4*(pn[1][j]+pn[2][j]+pn[1][j-1]+pn[1][j+1])
                +(pn[2][j]+pn[3][j]+pn[1][j-2]+pn[1][j+2]));
		}
        partvisc[1][MZ-1] = scale[MZ-1]*(12*pn[1][MZ-1]
			-4*(pn[1][MZ-1]+pn[2][MZ-1]+pn[1][MZ-2]+pn[1][MZ])
			+(pn[2][MZ-1]+pn[3][MZ-1]+pn[1][MZ-3]+pn[1][MZ]));
        partvisc[1][MZ] = scale[MZ]*(12*pn[1][MZ]
			-4*(pn[1][MZ]+pn[2][MZ]+pn[1][MZ-1]+pn[1][MZ])
			+(pn[2][MZ]+pn[3][MZ]+pn[1][MZ-2]+pn[1][MZ]));
// i = 2 terms																																	
        partvisc[2][1] = scale[1]*(12*pn[2][1]
			-4*(pn[1][1]+pn[3][1]+pn[2][1]+pn[2][2])
			+(pn[1][1]+pn[4][1]+pn[2][2]+pn[2][3]));
        partvisc[2][2] = scale[2]*(12*pn[2][2]
			-4*(pn[1][2]+pn[3][2]+pn[2][1]+pn[2][3])
			+(pn[1][2]+pn[4][2]+pn[2][1]+pn[2][4]));					
        for (j=3; j<=MZ-2; j++) { 
            partvisc[2][j] = scale[j]*(12*pn[2][j]
                    -4*(pn[1][j]+pn[3][j]+pn[2][j-1]+pn[2][j+1])
                    +(pn[1][j]+pn[4][j]+pn[2][j-2]+pn[2][j+2]));
        }
        partvisc[2][MZ-1] = scale[MZ-1]*(12*pn[2][MZ-1]
			-4*(pn[1][MZ-1]+pn[3][MZ-1]+pn[2][MZ-2]+pn[2][MZ])
			+(pn[1][MZ-1]+pn[4][MZ-1]+pn[2][MZ-3]+pn[2][MZ]));
		partvisc[2][MZ] = scale[MZ]*(12*pn[2][MZ]
			-4*(pn[1][MZ]+pn[3][MZ]+pn[2][MZ-1]+pn[2][MZ])
			+(pn[1][MZ]+pn[4][MZ]+pn[2][MZ-2]+pn[2][MZ]));
    }
// central grid
    for (i=ibegin; i<=iend; i++) {
		partvisc[i][1] = scale[1]*(12*pn[i][1]
				-4*(pn[i-1][1]+pn[i+1][1]+pn[i][1]+pn[i][2])
				+(pn[i-2][1]+pn[i+2][1]+pn[i][2]+pn[i][3]));
		partvisc[i][2] = scale[2]*(12*pn[i][2]
				-4*(pn[i-1][2]+pn[i+1][2]+pn[i][1]+pn[i][3])
				+(pn[i-2][2]+pn[i+2][2]+pn[i][1]+pn[i][4]));
        for (j=3; j<=MZ-2; j++) { 
            partvisc[i][j] = scale[j]*(12*pn[i][j]
				-4*(pn[i-1][j]+pn[i+1][j]+pn[i][j-1]+pn[i][j+1])
				+(pn[i-2][j]+pn[i+2][j]+pn[i][j-2]+pn[i][j+2]));
        }
		partvisc[i][MZ-1] = scale[MZ-1]*(12*pn[i][MZ-1]
				-4*(pn[i-1][MZ-1]+pn[i+1][MZ-1]+pn[i][MZ-2]+pn[i][MZ])
				+(pn[i-2][MZ-1]+pn[i+2][MZ-1]+pn[i][MZ-3]+pn[i][MZ]));
		partvisc[i][MZ] = scale[MZ]*(12*pn[i][MZ]
				-4*(pn[i-1][MZ]+pn[i+1][MZ]+pn[i][MZ-1]+pn[i][MZ])
				+(pn[i-2][MZ]+pn[i+2][MZ]+pn[i][MZ-2]+pn[i][MZ]));
    }	
	
// right hand terms MR-1, MR
	for (i=MR-1; i<=MR; i++) {
		partvisc[i][1] = scale[1]*(12*pn[i][1]
				-4*(pn[i-1][1]+pn[MR][1]+pn[i][1]+pn[i][2])
				+(pn[i-2][1]+pn[MR][1]+pn[i][2]+pn[i][3]));
		partvisc[i][2] = scale[2]*(12*pn[i][2]
				-4*(pn[i-1][2]+pn[MR][2]+pn[i][1]+pn[i][3])
				+(pn[i-2][2]+pn[MR][2]+pn[i][1]+pn[i][4]));
		for (j=3; j<=MZ-2; j++) { 
			partvisc[i][j] = scale[j]*(12*pn[i][j]
				-4*(pn[i-1][j]+pn[MR][j]+pn[i][j-1]+pn[i][j+1])
				+(pn[i-2][j]+pn[MR][j]+pn[i][j-2]+pn[i][j+2]));
	}
		partvisc[i][MZ-1] = scale[MZ-1]*(12*pn[i][MZ-1]
			-4*(pn[i-1][MZ-1]+pn[MR][MZ-1]+pn[i][MZ-2]+pn[i][MZ])
			+(pn[i-2][MZ-1]+pn[MR][MZ-1]+pn[i][MZ-3]+pn[i][MZ]));
		partvisc[i][MZ] = scale[MZ]*(12*pn[i][MZ]
			-4*(pn[i-1][MZ]+pn[MR][MZ]+pn[i][MZ-1]+pn[i][MZ])
			+(pn[i-2][MZ]+pn[MR][MZ]+pn[i][MZ-2]+pn[i][MZ]));
	
	}
	
	return;
}

/*######################################################################### 
/############################## artifiviscVX ############################*/
								 
void artifiviscVx(double **vxn, double **Vxartvisc, double *scale,  
		int iwindoleft, int iwindoright)
								 
// Vx is  symmetric about rigid surface and antisymmetric at r=0								 
{
	int i,j;
    int MR = (int) vxn[0][1]-1;
    int MZ = (int) vxn[0][0];
	//    printf("dimensions MR = %d MZ = %d \n",MR,MZ);
	
	int ibegin = int(fmax(3,iwindoleft));
	int iend = int(fmin(MR-1,iwindoright));

    if (iwindoleft<3) {
// i=2 terms (Vxartvisc[i=1][j] = 0)
        Vxartvisc[2][1] = scale[1]*(12*vxn[2][1]
			-4*(0+vxn[3][1]+vxn[2][1]+vxn[2][2])
			+(-vxn[2][1]+vxn[4][1]+vxn[2][2]+vxn[2][3]));
        Vxartvisc[2][2] = scale[2]*(12*vxn[2][2]
			-4*(0+vxn[3][2]+vxn[2][1]+vxn[2][3])
			+(-vxn[2][2]+vxn[4][2]+vxn[2][1]+vxn[2][4]));													 
        for (j=3; j<=MZ-2; j++) { 
			Vxartvisc[2][j] = scale[j]*(12*vxn[2][j]
				-4*(0+vxn[3][j]+vxn[2][j-1]+vxn[2][j+1])
				+(-vxn[2][j]+vxn[4][j]+vxn[2][j-2]+vxn[2][j+2]));
		}				
		Vxartvisc[2][MZ-1] = scale[MZ-1]*(12*vxn[2][MZ-1]
			-4*(0+vxn[3][MZ-1]+vxn[2][MZ-2]+vxn[2][MZ])
			+(-vxn[2][MZ-1]+vxn[4][MZ-1]+vxn[2][MZ-3]+vxn[2][MZ]));	
        Vxartvisc[2][MZ] = scale[MZ]*(12*vxn[2][MZ]
			-4*(0+vxn[3][MZ]+vxn[2][MZ-1]+vxn[2][MZ])
			+(-vxn[2][MZ]+vxn[4][MZ]+vxn[2][MZ-2]+vxn[2][MZ]));
    }
// central grid
    for (i=ibegin; i<=iend; i++) {
		Vxartvisc[i][1] = scale[1]*(12*vxn[i][1]
			-4*(vxn[i-1][1]+vxn[i+1][1]+vxn[i][1]+vxn[i][2])
			+(vxn[i-2][1]+vxn[i+2][1]+vxn[i][2]+vxn[i][3]));
		Vxartvisc[i][2] = scale[2]*(12*vxn[i][2]
			-4*(vxn[i-1][2]+vxn[i+1][2]+vxn[i][1]+vxn[i][3])
			+(vxn[i-2][2]+vxn[i+2][2]+vxn[i][1]+vxn[i][4]));													 
		for (j=3; j<=MZ-2; j++) { 
			Vxartvisc[i][j] = scale[j]*(12*vxn[i][j]
				-4*(vxn[i-1][j]+vxn[i+1][j]+vxn[i][j-1]+vxn[i][j+1])
				+(vxn[i-2][j]+vxn[i+2][j]+vxn[i][j-2]+vxn[i][j+2]));
		}												 
		Vxartvisc[i][MZ-1] = scale[MZ-1]*(12*vxn[i][MZ-1]
				-4*(vxn[i-1][MZ-1]+vxn[i+1][MZ-1]+vxn[i][MZ-2]+vxn[i][MZ])
				+(vxn[i-2][MZ-1]+vxn[i+2][MZ-1]+vxn[i][MZ-3]+vxn[i][MZ]));
												 
		Vxartvisc[i][MZ] = scale[MZ]*(12*vxn[i][MZ]
			-4*(vxn[i-1][MZ]+vxn[i+1][MZ]+vxn[i][MZ-1]+vxn[i][MZ])
			+(vxn[i-2][MZ]+vxn[i+2][MZ]+vxn[i][MZ-2]+vxn[i][MZ]));
	}																																			
// i= MR terms
	Vxartvisc[MR][1] = scale[1]*(12*vxn[MR][1]
			-4*(vxn[MR-1][1]+vxn[MR+1][1]+vxn[MR][1]+vxn[MR][2])
			+(vxn[MR-2][1]+vxn[MR+1][1]+vxn[MR][2]+vxn[MR][3]));
	Vxartvisc[MR][2] = scale[2]*(12*vxn[MR][2]
			-4*(vxn[MR-1][2]+vxn[MR+1][2]+vxn[MR][1]+vxn[MR][3])
			+(vxn[MR-2][2]+vxn[MR+1][2]+vxn[MR][1]+vxn[MR][4]));													 
	for (j=3; j<=MZ-2; j++) { 
		Vxartvisc[MR][j] = scale[j]*(12*vxn[MR][j]
			-4*(vxn[MR-1][j]+vxn[MR+1][j]+vxn[MR][j-1]+vxn[MR][j+1])
			+(vxn[MR-2][j]+vxn[MR+1][j]+vxn[MR][j-2]+vxn[MR][j+2]));
	}												 
	Vxartvisc[MR][MZ-1] = scale[MZ-1]*(12*vxn[MR][MZ-1]
			-4*(vxn[MR-1][MZ-1]+vxn[MR+1][MZ-1]+vxn[MR][MZ-2]+vxn[MR][MZ])
			+(vxn[MR-2][MZ-1]+vxn[MR+1][MZ-1]+vxn[MR][MZ-3]+vxn[MR][MZ]));
	
	Vxartvisc[MR][MZ] = scale[MZ]*(12*vxn[MR][MZ]
			-4*(vxn[MR-1][MZ]+vxn[MR+1][MZ]+vxn[MR][MZ-1]+vxn[MR][MZ])
			+(vxn[MR-2][MZ]+vxn[MR+1][MZ]+vxn[MR][MZ-2]+vxn[MR][MZ]));
	
	return;
}


/*######################################################################### 
/############################## artifiviscVX ############################*/

void artifiviscVy(double **vyn, double **Vyartvisc, double *scale,  
				  int iwindoleft, int iwindoright)

// Vy is  anti-symmetric about rigid surface and symmetric at r=0
// so Vyartvisc[i][j]==0
{
	int i,j;
    int MR = (int) vyn[0][1];
    int MZ = (int) vyn[0][0]-1;
	//    printf("dimensions MR = %d MZ = %d \n",MR,MZ);
	
	int ibegin = int(fmax(3,iwindoleft));
    int iend = int(fmin(MR-2,iwindoright));
	
//i=1
	Vyartvisc[1][2] = scale[2]*(12*vyn[1][2]
		-4*(vyn[1][2]+vyn[2][2]+0+vyn[1][3])
		+(vyn[2][2]+vyn[3][2]-vyn[1][2]+vyn[1][4]));	
	for (j=3; j<=MZ-1; j++) { 
		Vyartvisc[1][j] = scale[j]*(12*vyn[1][j]
			-4*(vyn[1][j]+vyn[2][j]+vyn[1][j-1]+vyn[1][j+1])
			+(vyn[2][j]+vyn[3][j]+vyn[1][j-2]+vyn[1][j+2]));		
	}
	Vyartvisc[1][MZ] = scale[MZ]*(12*vyn[1][MZ]
		-4*(vyn[1][MZ]+vyn[2][MZ]+vyn[1][MZ-1]+vyn[1][MZ+1])
		+(vyn[2][MZ]+vyn[3][MZ]+vyn[1][MZ-2]+vyn[1][MZ]));	
//i=2
	Vyartvisc[2][2] = scale[2]*(12*vyn[2][2]
		- 4*(vyn[1][2]+vyn[3][2]+0+vyn[2][3])
		+ (vyn[1][2]+vyn[4][2]-vyn[2][2]+vyn[2][4]));		
	for (j=3; j<=MZ-1; j++) { 
		Vyartvisc[2][j] = scale[j]*(12*vyn[2][j]
			- 4*(vyn[1][j]+vyn[3][j]+vyn[2][j-1]+vyn[2][j+1])
			+ (vyn[1][j]+vyn[4][j]+vyn[2][j-2]+vyn[2][j+2]));	
	}
	Vyartvisc[2][MZ] = scale[MZ]*(12*vyn[2][MZ]
		- 4*(vyn[1][MZ]+vyn[3][MZ]+vyn[2][MZ-1]+vyn[2][MZ+1])
		+ (vyn[1][MZ]+vyn[4][MZ]+vyn[2][MZ-2]+vyn[2][MZ+1]));
		
// central grid
    for (i=ibegin; i<=iend; i++) {
		Vyartvisc[i][2] = scale[2]*(12*vyn[i][2]
			-4*(vyn[i-1][2]+vyn[i+1][2]+0+vyn[i][3])
			+(vyn[i-2][2]+vyn[i+2][2]-vyn[i][2]+vyn[i][4]));	
		for (j=3; j<=MZ-1; j++) { 
			Vyartvisc[i][j] = scale[j]*(12*vyn[i][j]
				-4*(vyn[i-1][j]+vyn[i+1][j]+vyn[i][j-1]+vyn[i][j+1])
				+(vyn[i-2][j]+vyn[i+2][j]+vyn[i][j-2]+vyn[i][j+2]));	
		}
		Vyartvisc[i][MZ] = scale[MZ]*(12*vyn[i][MZ]
			-4*(vyn[i-1][MZ]+vyn[i+1][MZ]+vyn[i][MZ-1]+vyn[i][MZ+1])
			+(vyn[i-2][MZ]+vyn[i+2][MZ]+vyn[i][MZ-2]+vyn[i][MZ]));
	}																																			

// right hand side, i=MR-1, MR
	
	for (i=MR-1; i<=MR; i++) {
		Vyartvisc[i][2] = scale[2]*(12*vyn[i][2]
				-4*(vyn[i-1][2]+vyn[MR][2]+0+vyn[i][3])
				+(vyn[i-2][2]+vyn[MR][2]-vyn[i][2]+vyn[i][4]));	
		for (j=3; j<=MZ-1; j++) { 
			Vyartvisc[i][j] = scale[j]*(12*vyn[i][j]
				-4*(vyn[i-1][j]+vyn[MR][j]+vyn[i][j-1]+vyn[i][j+1])
				+(vyn[i-2][j]+vyn[MR][j]+vyn[i][j-2]+vyn[i][j+2]));	
		}
		Vyartvisc[i][MZ] = scale[MZ]*(12*vyn[i][MZ]
				-4*(vyn[i-1][MZ]+vyn[MR][MZ]+vyn[i][MZ-1]+vyn[i][MZ+1])
				+(vyn[i-2][MZ]+vyn[MR][MZ]+vyn[i][MZ-2]+vyn[i][MZ]));
	}			
	return;
}


/*######################################################################### 
/########### viscterm4 - like viscterm but 4th order accurate ###########*/

void viscterm4(double **vxvisc, double **vyvisc, double **vxn,double **vyn,
			double *viscx, double *viscy, double *viscylr, double *viscylz, 
			double *viscylrr, double *rvec, double *rvecx, int iwindoleft, 
			int iwindoright)

// this code makes use of symmetries at the left (source) boundary
// vx is an odd function about r=0: vx[-dr][j]=  = -vx[dr][j], vx[i][j]=0
// vy is an even function about r=0: vy[-dr][j] = vy[dr][j]
// BY METHOD OF IMAGES AT THE RIGID SURFACE, z=0
// vx is an even function & vy is an odd function 

// I have included cylindrical terms

{
    int i,j;
    int MR = (int) vyn[0][1];
    int MZ = (int) vxn[0][0];
	double a = 4.0/3, b=-1.0/12.0, c= -5.0, q = -1.0/6;
// coefficients for 1st derivative : (1/12,-2/3,0,2/3,-1/12)
// are multiplied by 2 to correct the 2*Dt not accounted for in viscylr
	
//	printf("iwindoleft= %d, iwindoright = %d\n",iwindoleft,iwindoright);
	
    if (iwindoleft==2) {        // source edge
		// Vx viscosity, i=2 terms  -  vxvisc[1][j] = 0
        vxvisc[2][1] = viscx[1]*(c*vxn[2][1]
				+ a*(0+vxn[3][1]+vxn[2][1]+vxn[2][2])
                + b*(-vxn[2][1]+vxn[4][1]+vxn[2][2]+vxn[2][3]))
				+ viscylr[1]*(a*(vxn[3][1]-0)+q*(vxn[4][1]+vxn[2][1])
				- vxn[2][1]*viscylrr[2])/rvecx[2];	
        vxvisc[2][2] = viscx[2]*(c*vxn[2][2]
				+ a*(0+vxn[3][2]+vxn[2][1]+vxn[2][3])
                + b*(-vxn[2][2]+vxn[4][2]+vxn[2][1]+vxn[2][4]))
				+ viscylr[2]*(a*(vxn[3][2]-0)+q*(vxn[4][2]+vxn[2][2])
				- vxn[2][2]*viscylrr[2])/rvecx[2];
        for (j=3; j<=MZ-2; j++) {	// i = 2 term
            vxvisc[2][j] = viscx[j]*(c*vxn[2][j]
				+ a*(0+vxn[3][j]+vxn[2][j-1]+vxn[2][j+1])
                + b*(-vxn[2][j]+vxn[4][j]+vxn[2][j-2]+vxn[2][j+2]))
				+ viscylr[j]*(a*(vxn[3][j]-0)+q*(vxn[4][j]+vxn[2][j])
				- vxn[2][j]*viscylrr[2])/rvecx[2];
        }
        for (j=MZ-1; j<=MZ; j++) {		// assume uniform vxn above MZ
			vxvisc[2][j] = viscx[j]*(c*vxn[2][j]
					+ a*(0+vxn[3][j]+vxn[2][j-1]+vxn[2][MZ])
                    + b*(-vxn[2][j]+vxn[4][j]+vxn[2][j-2]+vxn[2][MZ]))
					+ viscylr[j]*(a*(vxn[3][j]-0)+q*(vxn[4][j]+vxn[2][j]) 
					- vxn[2][j]*viscylrr[2])/rvecx[2];
		}
	// Vy, i=1 terms  , use symmetries      
        vyvisc[1][2] = viscy[2]*(c*vyn[1][2]		
                + a*(vyn[1][2]+vyn[2][2]+0+vyn[1][3])
                + b*(vyn[2][2]+vyn[3][2]-vyn[1][2]+vyn[1][4]))
                + viscylz[2]*(a*(vyn[2][2]-vyn[1][2]) 
                + q*(vyn[3][2]-vyn[2][2]))/rvec[1];
        for (j=3; j<=MZ-1; j++) { 
            vyvisc[1][j] = viscy[j]*(c*vyn[1][j]
                + a*(vyn[1][j]+vyn[2][j]+vyn[1][j-1]+vyn[1][j+1])
                + b*(vyn[2][j]+vyn[3][j]+vyn[1][j-2]+vyn[1][j+2]))
                + viscylz[j]*(a*(vyn[2][j]-vyn[1][j])
                + q*(vyn[3][j]-vyn[2][j]))/rvec[1]; 
        }
        vyvisc[1][MZ] = viscy[MZ]*(c*vyn[1][MZ]
                + a*(vyn[1][MZ]+vyn[2][MZ]+vyn[1][MZ-1]+vyn[1][MZ+1])
                + b*(vyn[2][MZ]+vyn[3][MZ]+vyn[1][MZ-2]+vyn[1][MZ+1]))
                + viscylz[MZ]*(a*(vyn[2][MZ]-vyn[1][MZ])
                + q*(vyn[3][MZ]-vyn[2][MZ]))/rvec[1]; 
		// Vy, i=2 terms         
        vyvisc[2][2] = viscy[2]*(c*vyn[2][2]
                + a*(vyn[1][2]+vyn[3][2]+0+vyn[2][3])
                + b*(vyn[1][2]+vyn[4][2]-vyn[2][2]+vyn[2][4]))
                + viscylz[2]*(a*(vyn[3][2]-vyn[1][2])
                + q*(vyn[4][2]-vyn[1][2]))/rvec[2]; 
        for (j=3; j<=MZ-1; j++) { 
            vyvisc[2][j] = viscy[j]*(c*vyn[2][j]
                + a*(vyn[1][j]+vyn[3][j]+vyn[2][j-1]+vyn[2][j+1])
                + b*(vyn[1][j]+vyn[4][j]+vyn[2][j-2]+vyn[2][j+2]))
                + viscylz[j]*(a*(vyn[3][j]-vyn[1][j])
                + q*(vyn[4][j]-vyn[1][j]))/rvec[2];             
        }
        vyvisc[2][MZ] = viscy[MZ]*(c*vyn[2][MZ]
                + a*(vyn[1][MZ]+vyn[3][MZ]+vyn[2][MZ-1]+vyn[2][MZ+1])
                + b*(vyn[1][MZ]+vyn[4][MZ]+vyn[2][MZ-2]+vyn[2][MZ+1]))
                + viscylz[MZ]*(a*(vyn[3][MZ]-vyn[1][MZ])
                + q*(vyn[4][MZ]-vyn[1][MZ]))/rvec[2];
    }
	
    if (iwindoright==MR) {            // far edge, use 2nd order accuracy
        vxvisc[MR][1] = viscx[1]*(vxn[MR-1][1]+vxn[MR+1][1]+vxn[MR][2]
					- 3*vxn[MR][1]) + viscylr[1]*(vxn[MR+1][1]-vxn[MR-1][1]
					- vxn[MR][1]*viscylrr[MR])/rvecx[MR];
        for (j=2; j<=MZ-1; j++) { 
            vxvisc[MR][j] = viscx[j]*(vxn[MR-1][j]+vxn[MR+1][j]
					+ vxn[MR][j-1]+vxn[MR][j+1]-4*vxn[MR][j])
					+ viscylr[j]*(vxn[MR+1][j]-vxn[MR-1][j]
					- vxn[MR][j]*viscylrr[MR])/rvecx[MR];
        }
        vxvisc[MR][MZ] = viscx[MZ]*(vxn[MR-1][MZ] + vxn[MR+1][MZ]
					+ vxn[MR][MZ-1] - 3*vxn[MR][MZ])
					+ viscylr[MZ]*(vxn[MR+1][MZ]-vxn[MR-1][MZ]
					- vxn[MR][MZ]*viscylrr[MR])/rvecx[MR];
        for (j=2; j<=MZ; j++) {
            vyvisc[MR][j] = viscy[j]*(vyn[MR][j-1] + vyn[MR][j+1]
						+ vyn[MR-1][j] - 3*vyn[MR][j]) + 2*viscylz[j]
						*(vyn[MR][j]-vyn[MR-1][j])/rvec[MR];
        }
    }
    
    int ibegin = int(fmax(3,iwindoleft));
    int iendx = int(fmin(MR-1,iwindoright));
    int iendy = int(fmin(MR-2,iwindoright));
    
	// main window - vxvisc terms		
    for (i=ibegin; i<=iendx; i++) {
        vxvisc[i][1] = viscx[1]*(c*vxn[i][1]
			+ a*(vxn[i-1][1]+vxn[i+1][1]+vxn[i][1]+vxn[i][2])
			+ b*(vxn[i-2][1]+vxn[i+2][1]+vxn[i][2]+vxn[i][3]))
			+ viscylr[1]*(a*(vxn[i+1][1]-vxn[i-1][1])
			+ q*(vxn[i+2][1]-vxn[i-2][j])-vxn[i][1]*viscylrr[i])/rvecx[i]; 
        vxvisc[i][2] = viscx[2]*(c*vxn[i][2]
			+ a*(vxn[i-1][2]+vxn[i+1][2]+vxn[i][1]+vxn[i][3])
			+ b*(vxn[i-2][2]+vxn[i+2][2]+vxn[i][1]+vxn[i][4]))
			+ viscylr[2]*(a*(vxn[i+1][2]-vxn[i-1][2])
			+ q*(vxn[i+2][2]-vxn[i-2][2])-vxn[i][2]*viscylrr[i])/rvecx[i];		
        for (j=3; j<=MZ-2; j++) { 
            vxvisc[i][j] = viscx[j]*(c*vxn[i][j]
				+ a*(vxn[i-1][j]+vxn[i+1][j]+vxn[i][j-1]+vxn[i][j+1])
				+ b*(vxn[i-2][j]+vxn[i+2][j]+vxn[i][j-2]+vxn[i][j+2]))
				+ viscylr[j]*(a*(vxn[i+1][j]-vxn[i-1][j])
				+ q*(vxn[i+2][j]-vxn[i-2][j])
				-vxn[i][j]*viscylrr[i])/rvecx[i];				
        }
        for (j=MZ-1; j<=MZ; j++) {		// assume uniform vxn above MZ
			vxvisc[i][j] = viscx[j]*(c*vxn[i][j]
				+ a*(vxn[i-1][j]+vxn[i+1][j]+vxn[i][j-1]+vxn[i][MZ])
				+ b*(vxn[i-2][j]+vxn[i+2][j]+vxn[i][j-2]+vxn[i][MZ]))
				+ viscylr[j]*(a*(vxn[i+1][j]-vxn[i-1][j])
				+ q*(vxn[i+2][j]-vxn[i-2][j])-vxn[i][j]*viscylrr[i])/rvecx[i];
		}
    }
	
	// main window - vyvisc terms
	for (i=ibegin; i<=iendy; i++) {
		vyvisc[i][2] = viscy[2]*(c*vyn[i][2]
                + a*(vyn[i-1][2]+vyn[i+1][2]+0+vyn[i][3])
                + b*(vyn[i-2][2]+vyn[i+2][2]-vyn[i][2]+vyn[i][4]))
                + viscylz[2]*(a*(vyn[i+1][2]-vyn[i-1][2])
                + q*(vyn[i+2][2]-vyn[i-2][2]))/rvec[i];
		for (j=3; j<=MZ-1; j++) { 
			vyvisc[i][j] = viscy[j]*(c*vyn[i][j]
					+ a*(vyn[i-1][j]+vyn[i+1][j]+vyn[i][j-1]+vyn[i][j+1])
					+ b*(vyn[i-2][j]+vyn[i+2][j]+vyn[i][j-2]+vyn[i][j+2]))
					+ viscylz[j]*(a*(vyn[i+1][j]-vyn[i-1][j])
					+ q*(vyn[i+2][j]-vyn[i-2][j]))/rvec[i]; 
		}
		vyvisc[i][MZ] = viscy[MZ]*(c*vyn[i][MZ]
				+ a*(vyn[i-1][MZ]+vyn[i+1][MZ]+vyn[i][MZ-1]+vyn[i][MZ+1])
				+ b*(vyn[i-2][MZ]+vyn[i+2][MZ]+vyn[i][MZ-2]+vyn[i][MZ+1]))
				+ viscylz[MZ]*(a*(vyn[i+1][MZ]-vyn[i-1][MZ])
				+ q*(vyn[i+2][MZ]-vyn[i-2][MZ]))/rvec[i];
	}
    
	return;

}

//######################################################################### 


