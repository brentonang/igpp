#include<iomanip> 
#include <cmath>
#include <ctime>           // needed for timestamp
#include <cstdio>          // needed for outputBinary_2d

# include "Array_aux.hpp"
# include "fdtdBoundaries_aux.hpp"

/*######################################################################### 
 
 FUNCTIONS SPECIFIC TO THE TOP AND RIGHT PML BOUNDARY CONDITIONS
 AS WELL AS THE RIGID BOUNDARY CONDITION WHERE THERE'S TOPOGRAPHY
 
########################################################################### 
/############################## rightBound ##############################*/
// do computations for the right PML boundary

void rightBound(int iebc, int NZ, int NR, double *dtdxorz,
		double **vxbcr, double **vxbcrn, double **vybcr, double **vybcrn,
		double **pxbcr, double **pxbcrn, double **pybcr, double **pybcrn,
		double **pybcb, double **pybcbn, double **pxbcbn, double **vybcbn,
		double **vx, double **vxn, double **pn,double *cavx2,double *cbvx2,
		double *cwrit, double **cavxbcr, double **cbvxbcr, double **dapybcb, 
		double **dbpybcb, double **dapxbcr, double **dbpxbcr, 
		double *dtdxrhocsq, double *pcylterm, double *rvec)
				
{    
    int i,j;
	int iefb = NR+iebc;
    
	swap2d(vxbcr,vxbcrn); swap2d(vybcr,vybcrn);
	swap2d(pxbcr,pxbcrn); swap2d(pybcr,pybcrn);
	
//  Update Vx, right boundary
	for (i=2; i<=iebc; i++)
		for (j=1; j<=NZ; j++)
			vxbcr[i][j]=cavxbcr[i][j]*vxbcr[i][j]-cbvxbcr[i][j]*
				(pxbcrn[i][j]+pybcrn[i][j]-pxbcrn[i-1][j]-pybcrn[i-1][j]); 
//	connect to Vx in main grid  	
	for (j=1; j<=NZ; j++)  // apply BC
		vx[NR+1][j] = cavx2[j]*vx[NR+1][j] - cbvx2[j]*
			(pxbcrn[1][j]+pybcrn[1][j]-pn[NR][j]);
	
//	Update Vy, right boundary       
	for (i=1; i<=iebc; i++)
		for (j=2; j<=NZ; j++)           //vybcr[i][1]=0 @ rigid bdy
			vybcr[i][j] = vybcr[i][j] - dtdxorz[j]*(pxbcrn[i][j]+
				pybcrn[i][j] - pxbcrn[i][j-1]-pybcrn[i][j-1]);
//	connect to top boundary (top right corner)	
	for (i=1; i<=iebc; i++)        // connect to back boundary
		vybcr[i][NZ+1]=vybcr[i][NZ+1]-dtdxorz[NZ]*(pxbcbn[NR+i][1]+
			pybcbn[NR+i][1]-pxbcrn[i][NZ]-pybcrn[i][NZ]);

/*.........................................................................
// for whatever reason, corrections for wind DO NOT decrease the reflection
// corrections for wind advection for particle velocities vxbcb
	for (i=2; i<=iebc; i++)
		for (j=1; j<=NZ; j++)
			vxbcr[i][j]=vxbcr[i][j]-cwrit[j]*(vxbcr[i+1][j]-vxbcr[i-1][j]);
	
	for (j=1; j<=NZ; j++)  // apply BC
		vx[NR+1][j] = vx[NR+1][j] - cwrit[j]*(vx[NR+1][j]-vx[NR][j])/2.0;
 
	for (i=2; i<=iebc-1; i++)
		for (j=2; j<=NZ; j++)           //vybcr[i][1]=0 @ rigid bdy
			vybcr[i][j] = vybcr[i][j]-cwrit[j]*(vybcr[i+1][j]-vybcr[i-1][j]);
	for (j=2; j<=NZ; j++)  {
		vybcr[1][j] = vybcr[1][j]-cwrit[j]*(vybcr[2][j]-vybcr[1][j])/2.0;
		vybcr[iebc][j] = vybcr[iebc][j]-cwrit[j]*(vybcr[iebc][j]-vybcr[iebc-1][j])/2.0;
	}
 
	for (i=2; i<=iebc-1; i++)        // connect to back boundary
		vybcr[i][NZ+1]=vybcr[i][NZ+1]-cwrit[NZ]*(vybcr[i+1][NZ+1]-vybcr[i-1][NZ+1]); 
	vybcr[1][NZ+1]=vybcr[1][NZ+1]-cwrit[NZ]*(vybcr[2][NZ+1]-vybcr[1][NZ+1])/2.0; 
	vybcr[iebc][NZ+1]=vybcr[iebc][NZ+1]-cwrit[NZ]*(vybcr[iebc][NZ+1]-vybcr[iebc-1][NZ+1])/2.0; 
*/
	
// Update Py, right boundary
	for (i=1; i<=iebc; i++)        // upper right corner			   
		pybcb[NR+i][1]=dapybcb[NR+i][1]*pybcb[NR+i][1]-
			dbpybcb[NR+i][1]*(vybcbn[NR+i][2]-vybcrn[i][NZ+1]);	
// Update Px, right boundary   
	for (i=2; i<=iebc; i++)
		for (j=1; j<=NZ; j++)
			pxbcr[i][j]=dapxbcr[i][j]*pxbcr[i][j]
				-dbpxbcr[i][j]*(vxbcrn[i+1][j]-vxbcrn[i][j])
				-pcylterm[j]*(vxbcrn[i+1][j]+vxbcrn[i][j])/rvec[NR];
// connect to Vx, main grid
	for (j=1; j<=NZ; j++)
		pxbcr[1][j]=dapxbcr[1][j]*pxbcr[1][j]
			-dbpxbcr[1][j]*(vxbcrn[2][j]-vxn[NR+1][j])
			-pcylterm[j]*(vxbcrn[2][j]+vxn[NR+1][j])/rvec[NR];	
// Update py, right boundary
	for (i=1; i<=iebc; i++)
		for (j=1; j<=NZ; j++)
			pybcr[i][j] = pybcr[i][j]-dtdxrhocsq[j]*
				(vybcrn[i][j+1]-vybcrn[i][j]);

/*.........................................................................
// for whatever reason, corrections for wind DO NOT decrease the reflection
// corrections for wind (skip for test B, include for C)
	for (i=1; i<=iebc-1; i++)        // upper right corner			   
		pybcb[NR+i][1] = pybcb[NR+i][1] - cwrit[NZ]*(pybcb[NR+i+1][1]-pybcb[NR+i-1][1]);
//	pybcb[iefb][1] = pybcb[iefb][1] - cwrit[NZ]*(pybcb[iefb][1]-pybcb[iefb-1][1])/2.0;
	
// Update Px for wind @ right boundary   
	for (i=2; i<=iebc-1; i++)
		for (j=1; j<=NZ; j++)
			pxbcr[i][j]=pxbcr[i][j]-cwrit[j]*(pxbcr[i+1][j]-pxbcr[i-1][j]);	
//	for (j=1; j<=NZ; j++)
//		pxbcr[iefb][j]=pxbcr[iefb][j]-cwrit[j]*(pxbcr[iefb][j]-pxbcr[iefb-1][j])/2.0;		
	
// connect to Vx, corrected for wind	
	for (j=1; j<=NZ; j++)
		pxbcr[1][j]=pxbcr[1][j]-cwrit[j]*(pxbcr[2][j]-pxbcr[1][j])/2.0;			
	
// Update py for wind @ right boundary
	for (i=2; i<=iebc-1; i++)
		for (j=1; j<=NZ; j++)
			pybcr[i][j] = pybcr[i][j] - cwrit[j]*(pybcr[i+1][j]-pybcr[i-1][j]);
	for (j=1; j<=NZ; j++){
		pybcr[1][j] = pybcr[1][j] - cwrit[j]*(pybcr[2][j]-pybcr[1][j])/2.0;
//		pybcr[iebc][j] = pybcr[iebc][j] - cwrit[j]*(pybcr[iebc][j]-pybcr[iebc-1][j])/2.0;	
	}
*/  
    return;
}


/*######################################################################### 
/################################ topBound ##############################*/
// do computations for the top PML boundary

void topBound(int iwindoright, int irb, int jebc, int NZ,  
			  double **vxbcb, double **vxbcbn, double **vybcb, double **vybcbn,
			  double **pxbcb, double **pxbcbn, double **pybcb, double **pybcbn,
			  double **vy, double **pn,double **vyn,double *cavy2,double *cbvy2,
			  double *cwtop,double **cavxbcb, double **cbvxbcb,double **cavybcb, 
			  double **cbvybcb, double **dapxbcb, double **dbpxbcb, 
			  double **dapybcb, double **dbpybcb, double *pcylbcb)

{    
    int i,j;
    
    swap2d(vxbcb,vxbcbn); swap2d(vybcb,vybcbn);
    swap2d(pxbcb,pxbcbn); swap2d(pybcb,pybcbn);
    
//  Update Vx, top boundary
    for (i=2; i<=irb; i++)       // vxbcb[1][j]=0 from symmetry
        for (j=1; j<=jebc; j++)
            vxbcb[i][j] = cavxbcb[i][j]*vxbcb[i][j]-cbvxbcb[i][j]*
                (pxbcbn[i][j]+pybcbn[i][j]-pxbcbn[i-1][j]-pybcbn[i-1][j]);
    
//  Update Vy, top boundary
    for (i=1; i<=irb; i++)
        for (j=2; j<=jebc-1; j++)
            vybcb[i][j] = cavybcb[i][j]*vybcb[i][j]-cbvybcb[i][j]*
                (pxbcbn[i][j]+pybcbn[i][j]-pxbcbn[i][j-1]-pybcbn[i][j-1]);
    
// connect Vy to main grid    
    for (i=1; i<=iwindoright; i++)   
        vy[i][NZ+1] = cavy2[i]*vy[i][NZ+1]-cbvy2[i]*
            (pxbcbn[i][1]+pybcbn[i][1]-pn[i][NZ]);

// corrections for wind advection for particle velocities vxbcb
    for (i=2; i<=irb-1; i++)       // vxbcb[1][j]=0 from symmetry
        for (j=1; j<=jebc; j++)
            vxbcb[i][j] = vxbcb[i][j]-cwtop[i]* (vxbcb[i+1][j]-vxbcb[i-1][j]);
	
// now do correction for wind at i= irb	
	for (j=1; j<=jebc; j++)
		vxbcb[irb][j] = vxbcb[irb][j]-cwtop[irb]* (vxbcb[irb][j]-vxbcb[irb-1][j])/2.0;	
	
// wind correction vybcb	
    for (i=2; i<=irb-1; i++)
        for (j=2; j<=jebc-1; j++)
            vybcb[i][j] = vybcb[i][j]-cwtop[i]* (vybcb[i+1][j]-vybcb[i-1][j]);	
	
// do correction for wind at i=1, irb
	for (j=2; j<=jebc-1; j++){
		vybcb[1][j] = vybcb[1][j]-cwtop[1]* (vybcb[2][j]-vybcb[1][j])/2.0;
		vybcb[irb][j] = vybcb[irb][j]-cwtop[irb]* (vybcb[irb][j]-vybcb[irb-1][j])/2.0;
	}	
// wind correction vy 	
	for (i=2; i<=iwindoright-1; i++)   {
        vy[i][NZ+1] = vy[i][NZ+1]-cwtop[i]* (vyn[i+1][NZ+1]-vyn[i-1][NZ+1]);
	}
	
// now do correction for wind at i=1, iwindoright 
	vy[1][NZ+1] = vy[1][NZ+1]-cwtop[1]* (vyn[2][NZ+1]-vyn[1][NZ+1])/2.0;
	vy[iwindoright][NZ+1] = vy[iwindoright][NZ+1]-cwtop[iwindoright]*
					(vyn[iwindoright][NZ+1]-vyn[iwindoright-1][NZ+1])/2.0;	
    
//  Update Px, top boundary
    for (i=1; i<=irb; i++) 
        for (j=1; j<=jebc; j++)	
            pxbcb[i][j]=dapxbcb[i][j]*pxbcb[i][j]
                -dbpxbcb[i][j]*(vxbcbn[i+1][j]-vxbcbn[i][j])
                -pcylbcb[i]*(vxbcbn[i+1][j]+vxbcbn[i][j]);
//  Update py, top boundary
    for (i=1; i <= irb; i++)
        for (j=2; j <= jebc; j++)
            pybcb[i][j]=dapybcb[i][j]*pybcb[i][j]-dbpybcb[i][j]*
                (vybcbn[i][j+1]-vybcbn[i][j]);
	
// wind correction 
    for (i=2; i<=irb-1; i++) 
        for (j=1; j<=jebc; j++)	
            pxbcb[i][j]=pxbcb[i][j] - cwtop[i]* (pxbcb[i+1][j]-pxbcb[i-1][j]);
	
    for (i=2; i <= irb-1; i++)
        for (j=2; j <= jebc; j++)
            pybcb[i][j]=pybcb[i][j] - cwtop[i]* (pybcb[i+1][j]-pybcb[i-1][j]);	
	
// wind corrections at 1=1,irb  
	pxbcb[1][1]=pxbcb[1][1] - cwtop[1]* (pxbcb[2][1]-pxbcb[1][1])/2.0;
	pxbcb[irb][1]=pxbcb[irb][1] - cwtop[irb]* (pxbcb[irb][1]-pxbcb[irb-1][1])/2.0;
	for (j=2; j <= jebc; j++) {
		pxbcb[1][j]=pxbcb[1][j] - cwtop[1]* (pxbcb[2][j]-pxbcb[1][j])/2.0;
		pxbcb[irb][j]=pxbcb[irb][j] - cwtop[irb]* (pxbcb[irb][j]-pxbcb[irb-1][j])/2.0;
		pybcb[1][j]=pybcb[1][j] - cwtop[1]* (pybcb[2][j]-pybcb[1][j])/2.0;
		pybcb[irb][j]=pybcb[irb][j] - cwtop[irb]* (pybcb[irb][j]-pybcb[irb-1][j])/2.0;
	}

    
// another connection to the main grid    
    for (i=1; i<=iwindoright; i++)       // connect to Vx in main grid
        pybcb[i][1]=dapybcb[i][1]*pybcb[i][1]-dbpybcb[i][1]*
            (vybcbn[i][2]-vyn[i][NZ+1]);
	
// correction for wind at pybcb[i][1]  (XXXX should try setting to iwindoright
    for (i=2; i<=iwindoright-1; i++)       // connect to Vx in main grid
        pybcb[i][1]=pybcb[i][1]- cwtop[i]* (pybcb[i+1][1]-pybcb[i-1][1]);	
	pybcb[1][1]=pybcb[1][1]- cwtop[1]* (pybcb[2][1]-pybcb[1][1])/2.0;
    
    return;
}

/*######################################################################### 
/################################ PMLsetup ##############################*/

// set up all the matrices of constants for the absorbing boundaries
// sound speeds and density vary with altitude. no adjustment for wind
void PMLsetup(int iebc, int NR,int NZ,double dt,double dx,double *cc,
		double *ww, double *rho,double *cavx2,double *cbvx2,double *cavy2,
		double *cbvy2, double *cwtop, double *cwrit, double **cavxbcb,double **cbvxbcb,double **cavybcb,
		double **cbvybcb,double **dapxbcb,double **dbpxbcb,double **dapybcb,
		double **dbpybcb,double **cavxbcr,double **cbvxbcr,double **dapxbcr,
		double **dbpxbcr,double *pcylbcb, double **pcylbcr)
{
//  setup & fill arrays for Perfectly Matched Layer (PML) boundaries
//  at top & right PML regions (bottom is rigid surface; left is r=0);
    
    int i,j,orderbc,ibbc,jbbc,ibfbc,jebc,iefbc;
    double reflmax,y1,y2,sigmay,sigmays,ca1,cb1,da1,db1,x1,x2;
    double delbc,sigmam,bcfactor, plambda, sigmax,sigmaxs;
    double *bcfactorv, *plambdav;
	double cwtbnd, cwrbnd;
    
    dt = 2*dt;
    jebc=iebc;             //thickness of right & back PML regions
    reflmax=0.0000001; orderbc=2;
    ibbc=iebc+1; jbbc=jebc+1; iefbc=NR+iebc; ibfbc=iefbc+1;
    char ofilename[80];
    
// PML fields: - derived ultimately from Hagness code
   
//     BACK region --------------------------------------------------------
    
    delbc=iebc*dx;
    sigmam=-log(reflmax)*rho[NZ]*cc[NZ]*(orderbc+1)/(2*delbc);
    bcfactor = sigmam/(dx*(pow(delbc,(double)orderbc))*(orderbc+1));
    plambda = 1/(cc[NZ]*cc[NZ]*rho[NZ]);
    
    printf("cc[NZ]= %lf, rho[NZ]= %14.12f in PMLsetup \n", cc[NZ],rho[NZ]);
	printf("ww[NZ]= %lf in PMLsetup \n", ww[NZ]);

	double damp = 0.0;		// see eqn 23. of JASA v124,pp 1430-1441 (2008)
    for (i=1; i<=iefbc; i++) {
        cavybcb[i][jbbc]=1.0-damp;
        cbvybcb[i][jbbc]=0.0;
        pcylbcb[i] = (1-damp/2)*dt/(plambda*(double(i)-0.5)*dx);
    }
    
    for (j=1; j<=jebc; j++) {
        y1=(j-0.5)*dx;
        y2=(j-1.5)*dx;
        sigmay=bcfactor*(pow(y1,(double)(orderbc+1))-pow(y2,(double)(orderbc+1)));
        ca1=exp(-sigmay*dt/rho[NZ]);
        cb1=(1-ca1)/(sigmay*dx);
        
//            printf("ca1 = %lf, cb1 = %lf \n", ca1,cb1);
        
        for (i=1; i<=iefbc; i++) {
            cavybcb[i][j]=ca1;
            cbvybcb[i][j]=cb1;
        }
    }
    
    sigmay = bcfactor*pow((0.5*dx),(double)(orderbc+1));;
    ca1=exp(-sigmay*dt/rho[NZ]);
    cb1=(1-ca1)/(sigmay*dx);
	cwtbnd = ww[NZ]*dt/(2*dx);
	
//    printf(" cwtbnd = %lf, ca1 = %lf, cb2 = %lf \n", cwtbnd, ca1,cb1);
    
    for (i=1; i<=NR; i++) { 
        cavy2[i]=ca1;
        cbvy2[i]=cb1;
	}
	for (i=1; i<=iefbc; i++) { 
		cwtop[i] = cwtbnd;
    }
	
//	printf("cavy2\n"); printvector(cavy2,2); printf(" \n"); 
//	printf("cbvy2\n"); printvector(cbvy2,2); printf(" \n"); 
	
    for (j=1; j<=jebc; j++) {
        y1=j*dx;
        y2=(j-1)*dx;
        sigmay=bcfactor*(pow(y1,(double)(orderbc+1))-pow(y2,(double)(orderbc+1)));
        sigmays=sigmay*(plambda/rho[NZ]);
        da1=exp(-sigmays*dt/plambda);
        db1=(1-da1)/(sigmays*dx);
        for (i=1; i<=iefbc; i++) {
            dapybcb[i][j]=da1;
            dbpybcb[i][j]=db1;
            cavxbcb[i][j]=1;
            cbvxbcb[i][j]=dt/(rho[NZ]*dx);
            dapxbcb[i][j]=(1.0-damp);
            dbpxbcb[i][j]=(1.0-damp/2)*dt/(plambda*dx);
        }
    }
    
// RIGHT REGION --------------------------------------------------------------
    
    bcfactorv=zeros(NZ); plambdav=zeros(NZ);
    for (j=1; j<=NZ; j++) {
        sigmam=-log(reflmax)*rho[j]*cc[j]*(orderbc+1)/(2*delbc);
        bcfactorv[j] = sigmam/(dx*(pow(delbc,(double)orderbc))*(orderbc+1));
        plambdav[j] = 1/(cc[j]*cc[j]*rho[j]);
    }
    
    for (j=1; j<=NZ; j++) {
        cavxbcr[ibbc][j]=1.0;
        cbvxbcr[ibbc][j]=0.0;
    }
    
    for (i=1; i<=iebc; i++) {
        x1=(i-0.5)*dx;
        x2=(i-1.5)*dx;
        for (j=1; j<=NZ; j++) {
            sigmax=bcfactorv[j]*(pow(x1,(double)(orderbc+1))-pow(x2,(double)(orderbc+1)));
            ca1=exp(-sigmax*dt/rho[j]);
            cb1=(1-ca1)/(sigmax*dx);
            cavxbcr[i][j]=ca1;
            cbvxbcr[i][j]=cb1;
        }
        for (j=1; j<=jebc; j++) {
            cavxbcb[i+NR][j]=ca1;
            cbvxbcb[i+NR][j]=cb1;
        }
    }
	/*
     sprintf(ofilename,"cavybcb%2.2d.bin",0);
     outputBinary_2d(ofilename,cavybcb,1,1.0,1.0,1);
     sprintf(ofilename,"cbvybcb%2.2d.bin",0);
     outputBinary_2d(ofilename,cbvybcb,1,1.0,1.0,1);
	 */
    for (j=1; j<=NZ; j++) {
        sigmax = bcfactorv[j]*pow((0.5*dx),(double)(orderbc+1));
        ca1 = exp(-sigmax*dt/rho[j]);
        cb1 = (1-ca1)/(sigmax*dx);
		cwrbnd = ww[j]*dt/(2*dx);
        cavx2[j]=ca1;
        cbvx2[j]=cb1;
		cwrit[j] = cwrbnd;
    }
	
//	printf("cwrit\n"); printvector(cwrit,6); printf(" \n");	exit(-1);
    
    for (i=1; i<=iebc; i++) {
        x1=i*dx;
        x2=(i-1)*dx;
        for (j=1; j<=NZ; j++) {
            sigmax=bcfactorv[j]*(pow(x1,(double)(orderbc+1))-pow(x2,(double)(orderbc+1)));
            sigmaxs=sigmax*plambdav[j]/rho[j];
            da1=exp(-sigmaxs*dt/plambdav[j]);
            db1=(1-da1)/(sigmaxs*dx);
            dapxbcr[i][j]=da1;
            dbpxbcr[i][j]=db1;
            pcylbcr[i][j]=da1*dt*cc[j]*cc[j]*rho[j]/((double(NR+i)-0.5)*dx);
        }
        for (j=1; j<=jebc; j++) {
            dapxbcb[i+NR][j]=da1;
            dbpxbcb[i+NR][j]=db1;
        }
    }
    
	/*
     sprintf(ofilename,"dapybcb%2.2d.bin",0);
     outputBinary_2d(ofilename,dapybcb,1,1.0,1.0,1);
     sprintf(ofilename,"dbpybcb%2.2d.bin",0);
     outputBinary_2d(ofilename,dbpybcb,1,1.0,1.0,1);     
	 */   
    delete[] bcfactorv; delete[] plambdav;
    return;
}

/*#########################################################################
/################ SET UP INDICES FOR A RIGID BOUNDARY ###################*/

    void topoSetup(double topo[], double zvec[], int jz[], int jx[])
{
    int i, nr = topo[0];
    double dum, z1 = zvec[1], dxx = zvec[2]-zvec[1];
    printf("nr = %d, z1 = %lf, dxx = %lf\n", nr,z1,dxx);
// first compute jztopo vector - where vz is set to zero in applyTopo
    for (i=1; i<=nr; i++) {
        jz[i] = 1+round((topo[i]-z1)/dxx);
    }

    jx[1]=jz[1]; jx[nr+1]=jz[nr];
// now compute jxtopo vector - vx=0 for j<jxtopo
    for (i=2; i<=nr; i++) {
        jx[i] = fmax(jz[i],jz[i-1]);
    }
            
    return;
}

/*#########################################################################
/################# APPLY INDICES FOR A RIGID BOUNDARY ###################*/

void applyTopo(int *jztopo, double **vz, int *jxtopo, double **vx)
{
    int i,j, nr = jztopo[0];

// first force vz to be zero at jztopo indices
    for (i=1; i<=nr; i++) {
        for (j=1; j<=jztopo[i]; j++) {
            vz[i][j] = 0.0;
        }
    }
    
// now force vx to be zero at vx (could also do p about here
    for (i=2; i<=nr; i++) {
        for (j=1; j<=jxtopo[i]-1; j++) {
            vx[i][j] = 0.0;
        }
    }

    return;
}

/*
#########################################################################
 /################ SET UP INDICES FOR A RIGID BOUNDARY ###################

int **topoSetup(double topo[], double zvec[], int jz[], int &nij)
{
    int i,j, jmin,jmax, icount, nr = topo[0];
    double dum, z1 = zvec[1];
    double dxx = zvec[2]-zvec[1];
    printf("nr = %d, z1 = %lf, dxx = %lf\n", nr,z1,dxx);
    // first compute jztopo vector - where vz is set to zero in applyTopo
    for (j=1; j<=nr; j++) {
        jz[j] = 1+round((topo[j]-z1)/dxx);
    }
    // now find nij - the length of the 2 x nij output array
    nij = 0;
    for (j=2; j<=nr; j++) {
        nij = nij + abs(jz[j]-jz[j-1]);
    }
    //    printf("nij = %d\n", nij);
    // now set ijpair array - where vx is set to zero in applyTopo
    int **ijpair = izeros(2,nij);
    icount = 0;
    for (i=2; i<=nr; i++) {
        if (jz[i] != jz[i-1]) {
            jmin = fmin(jz[i],jz[i-1]);
            jmax = fmax(jz[i],jz[i-1])-1;
            for (j=jmin; j<=jmax; j++) {
                icount = icount+1;
                ijpair[1][icount] = i;
                ijpair[2][icount] = j;
            }
        }
    }
    if (icount !=nij){
        printf(" FATAL ERROR \n");
        printf(" problem in topoSetup \n");
        exit(-1);
    }
    
    return ijpair;
}

#########################################################################
 /################# APPLY INDICES FOR A RIGID BOUNDARY ###################

void applyTopo(int *jztopo, double **vz,int **ijpairs,double **vx,int nij)
{
    int i,j, k, nr = jztopo[0];
    
    // first force vz to be zero at jztopo indices
    for (i=1; i<=nr; i++) {
        j = jztopo[i];
        vz[i][j] = 0.0;
    }
    
    // now force vx to be zero at ijpairs
    for (k=1; k<=nij; k++) {
        i = ijpairs[1][k];
        j = ijpairs[2][k];
        vx[i][j] = 0.0;
    }
    
    return;
}

*/
