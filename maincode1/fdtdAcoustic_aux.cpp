#include <cmath>
#include <cstdio>
#include <ctime>           // needed for timestamp          
#include <iomanip> 

#include "Array_aux.hpp"
#include "fdtdWave_aux.hpp"
#include "fdtdAcoustic_aux.hpp"

/*######################################################################### 
 
 FUNCTIONS SPECIFIC TO FDTD Acoustic Propagation
 
###############################d#########################################*/

/*######################################################################### 
/### acoustSetup - set constants for time step loop - 1D env model ###*/

// there's a trade off between the memory required for the constants & the
// computational speed (more multiplication/division operations slows it down)
// SO, in this version, 1/r term is taken outside of cylindrical terms so that
// 1D vectors can be used instead of 2D arrays for pcylterm and visccylatt
//
// dt/(dx*rho) have to be computed separately for vx and vz since they're
// very slightly offset in altitude


// XXX; still need to add wind shear


void acoustSetup(double *Vxacterm,double *Vyacterm, double *Pacterm,
			double *pcylterm, double *advecterm, double *shearterm, 
			double *rho, double *c1d, double *w1d, double dt, double dx) {
   int i,j;
	int NZ = (int) rho[0];
    
   for (j = 1; j <= NZ; j++) {
      Vxacterm[j] = 2*dt/(dx*rho[j]);            
      Pacterm[j] = 2*dt*c1d[j]*c1d[j]*rho[j]/dx;  
      pcylterm[j] = dt*c1d[j]*c1d[j]*rho[j];    // p cyl term-needs 1/r
		advecterm[j] = dt*w1d[j]/dx; 
   }
	Vyacterm[1] = Vxacterm[1]; Vyacterm[NZ] = Vxacterm[NZ];
	for (j = 2; j <= NZ - 1; j++) {
		Vyacterm[j] = 2*dt/(dx*(rho[j-1]+rho[j])/2);
		shearterm[j] = dt * (w1d[j+1]-w1d[j-1])/(4*dx);	// average over 4 adjacent vz nodes
	}
	shearterm[1] = shearterm[2]; shearterm[NZ] = shearterm[NZ-1];
   return; // return not necessary
}

/*######################################################################### 
/###### acoust_setupN - constants needed for BC in nonlinear case #######*/

void acoust_setupN(double *dtdxorz, double *dtdxrhocsq, double *pcylterm, 
				double *rho, double *c1d, double dt,double dx) {
	int j;
	int NZ = (int) rho[0];
	
   for (j = 1; j <= NZ; j++) {
      dtdxrhocsq[j] = 2*dt*c1d[j]*c1d[j]*rho[j]/dx;  //P term for BC
      pcylterm[j] = dt*c1d[j]*c1d[j]*rho[j];         //P cyl term for BC
   }
	dtdxorz[1] = 2*dt/(dx*rho[1]); dtdxorz[NZ] = 2*dt/(dx*rho[NZ]);
	for (j = 2; j <= NZ-1; j++) dtdxorz[j] = 2*dt/(dx*(rho[j-1]+rho[j])/2);
	return;
}

/*######################################################################### 
/#################################  TstepOne ############################*/

void TstepOne(double **vx, double **vy, double **p, double **vxn, 
			double **vyn, double **pn, double *Vxacterm, double *Vyacterm, 
			double *Pacterm, double *pcylterm, double *rvec, int *jztopo,
            int *jxtopo, double dt, double dx, int widx) {
   // allows for topography
   int i,j, j1;
	int NZ = (int) pn[0][0];
	
   // do forward steps for first time step; Note time step is over 1*dt
   // V fields will be separated in time from P by 1/2 time steps
   // V fields at time = dt/2, start at i=2 for vx[1][j]=0
   for (i = 2; i <= widx; i++) {
      j1 = jxtopo[i];    //start at j=j1, not j=1
      for (j = j1; j <= NZ; j++) {
         vx[i][j] = vxn[i][j]-Vxacterm[j]*(pn[i][j]-pn[i-1][j])/2;
      }   
   }
    
   for (i = 1; i <= widx; i++) {
      // rigid surface v[i][jztopo[i]]=0, start @ j=j1, not j=2 (vy[i][NZ]=0 @ source)
      j1 = jztopo[i]+1;
      for (j = j1; j <= NZ; j++) {
         vy[i][j] = vyn[i][j]-Vyacterm[j]*(pn[i][j]-pn[i][j-1])/2;
      }   
		
      // P field at time dt,      start at j=j1, not j=1
      j1 = jztopo[i];
      for (j = j1; j <= NZ; j++) {
         p[i][j] = pn[i][j] - (Pacterm[j]/2)*(vy[i][j+1]-vy[i][j]+ 
			   vx[i+1][j] -vx[i][j]) - (pcylterm[j]/2)*
			   (vx[i+1][j]+vx[i][j])/rvec[i];
      }	
   }
   return;
}

/*######################################################################### 
/####################### Tstep_OneN nonlinear ###########################*/

void Tstep_OneN(double **vx,double **vy,double **r,double **vxn,
				double **vyn, double **rn, double **pn, double *rho,  
				double *cc, double dt, double dx, int widx) {
   /* similar to TstepOne except that solution is for r and v, not p and v
    include "self-advection" terms (these are ==0 for v for the 1st time-steps)
    As for linear Tstep_first, I neglect attenuation terms
   */
   int i,j;
	int NZ = (int) r[0][0];
	double dtdx = dt/dx;
	double dtd2x = dtdx/2, dtd4x = dtdx/4;
   for (i = 2; i <= widx; i++) {           // start at 2 for vx[1][j]=0
      for (j = 1; j <= NZ; j++) {
         vx[i][j] = vxn[i][j] - dtdx*
		 	  (pn[i][j]-pn[i-1][j])/(rho[j]+(rn[i][j]+rn[i-1][j])/2.);
      }   // without any vx[NR][j] values, this will be a rigid boundary
   } 
   for (i = 1; i <= widx; i++) {
		// rigid surface v[i][1] = 0, (also vy[i][NZ+1]=0, probably not a problem)
      for (j=2; j<=NZ; j++) {
         vy[i][j] = vyn[i][j] - dtdx*
		    	(pn[i][j]-pn[i][j-1])/(rho[j]+(rn[i][j]+rn[i][j-1])/2.);
      }
	}
	for (i = 1; i <= widx; i++) {
		for (j = 1; j <= NZ; j++) {		
         r[i][j] = rn[i][j]-(rho[j]+rn[i][j])*
				(dtdx*(vy[i][j+1]-vy[i][j]+ vx[i+1][j]-vx[i][j])
				 +dt*(vx[i+1][j]+vx[i][j])/(2*((double)i-0.5)*dx));	// cyl term 
		}
	}
	// add self-advection terms to the density matrix 	
	for (i = 2; i <= widx; i++) {
		r[i][1] = r[i][1] - dtd4x*((vx[i+1][1]+vx[i][1])*(rn[i+1][1]-rn[i-1][1]))
			-  dtd2x* ( (vy[i][1]+vy[i][2]) * (rn[i][2]-rn[i][1]) );
		for (j = 2; j <= NZ-1; j++) {		
			r[i][j] = r[i][j] - dtd4x*( (vx[i+1][j]+vx[i][j])*(rn[i+1][j]-rn[i-1][j])		
					+ (vy[i][j]+vy[i][j+1]) * (rn[i][j+1]-rn[i][j-1]) );	
      }
		r[i][NZ] = r[i][NZ] - dtd4x*((vx[i+1][NZ]+vx[i][NZ])*(rn[i+1][NZ]-rn[i-1][NZ]))		
				-dtd2x* ( (vy[i][NZ]+vy[i][NZ+1]) * (rn[i][NZ]-rn[i][NZ-1]) );
   }
	// now do self-advection terms at i==1	
	r[1][1] = r[1][1] - dtd2x*( vx[2][1]*(rn[2][1]-rn[1][1]) 
							   + (vy[1][1]+vy[1][2])*(rn[1][2]-rn[1][1]));	
	for (j = 2; j <= NZ-1; j++) {
		r[1][j] = r[1][j] - dtd2x*( vx[2][1]*(rn[2][j]-rn[1][j]) )		
				- dtd4x*(vy[1][j]+vy[1][j+1])*(rn[1][j+1]-rn[1][j-1]);	
	}
	r[1][NZ] = r[1][NZ] - dtd2x*( vx[2][NZ]*(rn[2][NZ]-rn[1][NZ]) 
				+ (vy[1][NZ]+vy[1][NZ+1])*(rn[1][NZ]-rn[1][NZ-1]));			
   return;
}

/*######################################################################### 
/#################### acoustic terms - 2nd order accuracy ###############*/

void acoustTerm(double **pacoust, double **vxacoust, double **vyacoust, 
			double **pn, double **vxn, double **vyn, double *Vxacterm,
			double *Vyacterm, double *Pacterm, double *pcylterm, 
			double *wterm, double *swterm, double *rvec, int iwindoleft, 
			int iwindoright, int jwindotop, bool ileft, bool iright) {
   // would it be faster to have the outer loop over j? i.e. row order vs. column order?
   int i,j, ir1, im1 = iwindoleft-1; 
   int MR = (int) pn[0][1];
   int NZ = (int) pn[0][0]; 
	double damp = 0.01;
	
	if (iright) {	// if computation window HAS reached right end
		ir1 = iwindoright-1;
		i = iwindoright;	// deal with i+1 indices at right for wind
		
		vxacoust[i][1] = -Vxacterm[1]*(pn[i][1]-pn[i-1][1]) 
						- wterm[1]*(vxn[i+1][1]-vxn[i-1][1])
						- swterm[1]*(vyn[i][2]+vyn[i-1][2]);
		pacoust[i][1] = -Pacterm[1]*(vyn[i][2]+vxn[i+1][1]-vxn[i][1])
						- pcylterm[1]*(vxn[i+1][1]+vxn[i][1])/rvec[i]
						- wterm[1]*(pn[i][1]-pn[i-1][1])/2.0;
		
		for (j = 2; j <= jwindotop; j++) {			// main grid
			vxacoust[i][j] = -Vxacterm[j]*(pn[i][j]-pn[i-1][j])
				- wterm[j]*(vxn[i+1][j]-vxn[i-1][j]) - swterm[j]* 
				(vyn[i][j+1] + vyn[i][j] + vyn[i-1][j+1] + vyn[i-1][j]); 
			vyacoust[i][j] = -Vyacterm[j]*(pn[i][j]-pn[i][j-1])
							- wterm[j]*(vyn[i][j]-vyn[i-1][j])/2.0;
			pacoust[i][j] = -Pacterm[j]*(vyn[i][j+1]-vyn[i][j]+ 
							vxn[i+1][j] -vxn[i][j]) - pcylterm[j]*
							(vxn[i+1][j]+vxn[i][j])/rvec[i]
							- wterm[j]*(pn[i][j]-pn[i-1][j])/2.0;	
		}		
	} else {		// if computation window has NOT reached right end
		ir1 = iwindoright;
	}
	for (i = iwindoleft; i <= ir1; i++) {	// terms at Earth bdy	
		vxacoust[i][1] = -Vxacterm[1]*(pn[i][1]-pn[i-1][1]) 
						- wterm[1]*(vxn[i+1][1]-vxn[i-1][1])   
						- swterm[1]*(vyn[i][2] + vyn[i-1][2]);
		pacoust[i][1] = -Pacterm[1]*(vyn[i][2]+vxn[i+1][1]-vxn[i][1])
				- pcylterm[1]*(vxn[i+1][1]+vxn[i][1])/rvec[i]
				- wterm[1]*(pn[i+1][1]-pn[i-1][1]);			
		for (j = 2; j <= jwindotop; j++) {			// main grid
			vxacoust[i][j] = -Vxacterm[j]*(pn[i][j]-pn[i-1][j])
				- wterm[j]*(vxn[i+1][j]-vxn[i-1][j])- swterm[j]*   
				(vyn[i][j+1] + vyn[i][j] + vyn[i-1][j+1] + vyn[i-1][j]);
			vyacoust[i][j] = -Vyacterm[j]*(pn[i][j]-pn[i][j-1])
						- wterm[j]*(vyn[i+1][j]-vyn[i-1][j]);
			pacoust[i][j] = -Pacterm[j]*(vyn[i][j+1]-vyn[i][j]+ 
						vxn[i+1][j] -vxn[i][j]) - pcylterm[j]*
						(vxn[i+1][j]+vxn[i][j])/rvec[i]
						- wterm[j]*(pn[i+1][j]-pn[i-1][j]);	
		}
	}
   if (ileft) {  // source region
      for (j = 1; j <= jwindotop; j++) {             
         pacoust[1][j] = -Pacterm[j]*(vyn[1][j+1]-vyn[1][j]  
                       + vxn[2][j])-pcylterm[j]*vxn[2][j]/rvec[1]
		    				  - wterm[j]*(pn[2][j]-pn[1][j])/2.0;
      }
      for (j=2; j<=jwindotop; j++) {
         vyacoust[1][j] = -Vyacterm[j]*(pn[1][j]-pn[1][j-1])
							   - wterm[j]*(vyn[2][j]-vyn[1][j])/2.0;
      }
   } else {	 // the signal has left the source region
      vxacoust[im1][1] = 0.98*vxacoust[iwindoleft][1];
      pacoust[im1][1] = -Pacterm[1]*(vyn[im1][2]+vxn[im1+1][1]-
				vxn[im1][1])- pcylterm[1]*(vxn[im1+1][1]+vxn[im1][1])/rvec[im1]
				- wterm[1]*(pn[iwindoleft][1]-pn[im1][1])/2.0;
      for (j = 2; j <= jwindotop; j++) {			// main grid
         vxacoust[im1][j] = 0.98*vxacoust[iwindoleft][j];		
         vyacoust[im1][j] = -Vyacterm[j]*(pn[im1][j]-pn[im1][j-1])
							- wterm[j]*(vyn[iwindoleft][j]-vyn[im1][j])/2.0;
         pacoust[im1][j] = -Pacterm[j]*(vyn[im1][j+1]-vyn[im1][j]+ 
							vxn[im1+1][j] -vxn[im1][j]) - pcylterm[j]*
							(vxn[im1+1][j]+vxn[im1][j])/rvec[im1]
							- wterm[j]*(pn[iwindoleft][j]-pn[im1][j])/2.0;	
      }		
	}    
   return;
}

/*######################################################################### 
/#################### acoustic terms - 2nd order accuracy ###############*/

void acoustTerm4(double **pacoust, double **vxacoust, double **vyacoust, 
                double **pn, double **vxn, double **vyn, double *Vxacterm,
                double *Vyacterm, double *Pacterm, double *pcylterm,
                double *rvec, int iwindoleft, int iwindoright, bool ileft) {

   // 4th order accuracy for a regular staggered grid, (from Fornberg)
   // df/dx = 9/8*(f(x+h/2)-f(x-h/2))-1/12*(f(x+3h/2)-f(x-3h/2))

   // this code makes use of symmetries at the left (source) boundary
   // vx is an odd function about r=0: vx[-dr][j]=  = -vx[dr][j], vx[i][j]=0
   // vy and p are even functions about r=0: vy[-dr][j] = vy[dr][j]
   // BY METHOD OF IMAGES AT THE RIGID SURFACE, z=0
   // vx & p are even functions & vy is an odd function 

   // cylindrical term is still 2nd order?
   int i,j;
   int MR = (int) pn[0][1];
   int NZ = (int) pn[0][0];
   double qa = 9/8, qb = -1/12;

   // main grid
   int ibegin = int(fmax(3,iwindoleft));
   int iendx = int(fmin(MR-1,iwindoright));
	int ibeginy = ileft? 1:iwindoleft;	
        
   for (i = ibegin; i <= iendx; i++) {	
      for (j = 1; j <= NZ; j++) {
         vxacoust[i][j] = -Vxacterm[j]*(qa*(pn[i][j]-pn[i-1][j])
                                + qb*(pn[i+1][j]-pn[i-2][j]));
      }
   }

   // all Vy terms
   for (i = ibeginy; i <= iwindoright; i++) {
      // j = 2 term, use mirror image symmetry at rigid surface
      vyacoust[i][2] = -Vyacterm[2]*(qa*(pn[i][2]-pn[i][1])
                                + qb*(pn[i][3]-pn[i][1]));
      for (j = 3; j <= NZ-1; j++) {
         vyacoust[i][j] = -Vyacterm[j]*(qa*(pn[i][j]-pn[i][j-1])
                                + qb*(pn[i][j+1]-pn[i][j-2]));
      }
      // j = NZ term, 2nd order accuracy
      vyacoust[i][NZ] = -Vyacterm[NZ]*(pn[i][NZ]-pn[i][NZ-1]);
	}
	
   // p terms
	for (i = iwindoleft; i <= iendx; i++) {
      // j=1, use symmetry at rigid surface
      pacoust[i][1] = -Pacterm[1]*(qa*(vyn[i][2]-0+ 
                    vxn[i+1][1]-vxn[i][1]) + qb*(vyn[i][3]+vyn[i][2]+ 
                    vxn[i+2][1]-vxn[i-1][1])) - pcylterm[1]*
                    (vxn[i+1][1]+vxn[i][1])/rvec[i];
      for (j = 2; j <= NZ-1; j++) {
         pacoust[i][j] = -Pacterm[j]*(qa*(vyn[i][j+1]-vyn[i][j]+ 
                 vxn[i+1][j]-vxn[i][j]) + qb*(vyn[i][j+2]-vyn[i][j-1]+ 
                 vxn[i+2][j]-vxn[i-1][j])) - pcylterm[j]*
                 (vxn[i+1][j]+vxn[i][j])/rvec[i];
      }
      // j=NZ, 2nd order accuracy
      pacoust[i][NZ] = -Pacterm[NZ]*(vyn[i][NZ+1]-vyn[i][NZ]+ 
               vxn[i+1][NZ] -vxn[i][NZ]) - pcylterm[NZ]*
               (vxn[i+1][NZ]+vxn[i][NZ])/rvec[i];	
   }                                        
   if (ileft) {	// source region
		for (j = 1; j <= NZ; j++) {		//pn is symmetric about r=0
         vxacoust[2][j] = -Vxacterm[j]*(qa*(pn[2][j]-pn[1][j])
							+ qb*(pn[3][j]-pn[1][j]));
		}
      // p terms 
      pacoust[1][1] = -Pacterm[1]*(qa*(vyn[i][2]-0+ vxn[2][1]-0) 
                    + qb*(vyn[1][3]+vyn[1][2] + vxn[3][1]+vxn[2][1])) 
                    - pcylterm[1]* (vxn[2][1]+0)/rvec[1];
      for (j = 2; j <= NZ-1; j++) {
         pacoust[1][j] = -Pacterm[j]*(qa*(vyn[1][j+1]-vyn[1][j] 
                 + vxn[2][j]-0) + qb*(vyn[1][j+2]-vyn[1][j-1] 
                 + vxn[3][j]+vxn[2][j])) - pcylterm[j]*
                 (vxn[2][j]+0)/rvec[1];   
      }
      // j=NZ, 2nd order accuracy
      pacoust[1][NZ] = -Pacterm[NZ]*(vyn[1][NZ+1]-vyn[1][NZ]+ 
                  vxn[1][NZ]-0) - pcylterm[NZ]*(vxn[2][NZ]+0)/rvec[i];        
   }	
    
   if (iwindoright == MR) {		// 2nd order accuracy for vx and p
      for (j=1; j<=NZ; j++) {         
         vxacoust[MR][j] = -Vxacterm[j]*(pn[MR][j]-pn[MR-1][j]);
			pacoust[MR][j] = -Pacterm[j]*(vyn[MR][j+1]-vyn[MR][j]+ 
								vxn[MR+1][j] -vxn[MR][j]) - pcylterm[j]*
								(vxn[MR+1][j]+vxn[MR][j])/rvec[i];	
      }		    
   }    
   return;
}

/*######################################################################### 
/#################### acoustic terms - 2nd order accuracy ###############*/

void NonlinAcoust(double **rnacoustnl, double **vxacoustnl, 
			double **vyacoustnl, double **pn, double **rn, double **vxn, 
			double **vyn, double *rho, double *rvec,  double d2tdx, 
			double d4tdx, double d2t, double dx, int iwindoleft, 
			int iwindoright, bool ileft) {
   int i,j;
   int MR = (int) rn[0][1];
   int NZ = (int) rn[0][0];
    
	for (i = iwindoleft; i <= iwindoright; i++) {
		vxacoustnl[i][1] = -d2tdx*(pn[i][1]-pn[i-1][1])/
				(rho[1]+(rn[i][1]+rn[i-1][1])/2.);		
		rnacoustnl[i][1] = -d2t*(rho[1]+rn[i][1])* 
				((vyn[i][2]-vyn[i][1]+vxn[i+1][1]-vxn[i][1])/dx 
				 +0.5*(vxn[i+1][1]+vxn[i][1])/rvec[i]);
		for (j = 2; j <= NZ; j++) {			// main grid
			vxacoustnl[i][j] = -d2tdx*(pn[i][j]-pn[i-1][j])/
				(rho[j]+(rn[i][j]+rn[i-1][j])/2.);	

			vyacoustnl[i][j] = -d4tdx*(pn[i][j]-pn[i][j-1])/
				(rho[j]+rho[j-1]+rn[i][j]+rn[i][j-1]);
			
			rnacoustnl[i][j] = -d2t*(rho[j]+rn[i][j])*
				((vyn[i][j+1]-vyn[i][j]+vxn[i+1][j]-vxn[i][j])/dx 
				 +0.5*(vxn[i+1][j]+vxn[i][j])/rvec[i]);					
		}			
	}
	if (ileft) {						// source region	
		for (j = 1; j <= NZ; j++) {             // 1st diff for wind term
			rnacoustnl[1][j] = -d2t*(rho[j]+rn[1][j])*
				((vyn[1][j+1]-vyn[1][j] + vxn[2][j])/dx 
				 +0.5*vxn[2][j]/rvec[1]);    //vxn[i][j]=0 @ boundary   
		}
		for (j = 2; j <= NZ; j++) {
			vyacoustnl[1][j] = -d4tdx*(pn[1][j]-pn[1][j-1])/
				(rho[j]+rho[j-1]+rn[1][j]+rn[1][j-1]);
		}
	}	
   return;
}
//######################################################################### 