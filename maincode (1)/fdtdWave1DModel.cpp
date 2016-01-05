  #include <iomanip>
#include <cmath>
#include <ctime>           // needed for timestamp
#include <cstdio>          // needed for outputBinary_matrix

# include "fdtdWave_aux.hpp"
# include "fdtdBoundaries_aux.hpp"
# include "fdtdViscosity_aux.hpp"
# include "fdtdAcoustic_aux.hpp"
# include "Array_aux.hpp"
# include "IO_utils.hpp"

/*#########################################################################
/##################### more function prototypes   #######################*/

const double grav = 9.8;

void acoustic_1Dmodel(int rmax,int zmin,int zmax,float fsrc,int srcz,
					  float Tend);

/*#########################################################################
/#############################  MAIN   ##################################*/

int main()
{    
// SET THE FOLLOWING CONSTANTS BEFORE COMPILATION
    
    const int zmax = 19000;			// maximum altitude in meters
    const int zmin = 1260;			// altitude of rigid surface in meters
    const int srcalt = 1261;		// source altitude in meters
    const int rmax = 32000;		// maximum range in meters
    const float fsrc_max = 5.0;     // maximum source frequency
	const float Tfin = 88.001;		// solve to time Tfin (seconds)
    
// DIMENSION = the dimension of the environmental model
//          (sound & wind speeds, density & pressure)
//      dimension=1: the environmental model varies in altitude
//      dimension=2: the environmental model varies in altitude & range
    
    const int dimension = 1;
    
    if (srcalt<zmin || zmax<=srcalt) {         
        printf(" FATAL ERROR \n");
        printf(" max alt=%d, srcalt=%d, min alt = %d\n",zmax,srcalt,zmin);
		printf(" SUGGESTION: fix s.t. zmin <= srcalt < zmax & recompile \n");
        exit(-1);  
    }  
//.........................................................................
    
    switch (dimension) {
        case 1:
            printf("\n 3D Acoustic propagation in a 1D model \n");
            printf(" sound & wind speeds & density vary with altitude\n");
            acoustic_1Dmodel(rmax,zmin,zmax,fsrc_max,srcalt,Tfin);
            break;
        case 2:
//            acoustic_2Dmodel(rmax,zmin,zmax,fsrc_max,srcalt,Tfin);
            printf(" Model dimension=2  DIY \n");
            break;
        default:
            printf("The model dimension must be either 1, or 2\n");
            exit(1);            
    }    
    return 0;
}

/*#########################################################################
/############################ acoustic_1D model #########################*/
void acoustic_1Dmodel(int rmax,int zmin,int zmax,float fsrc,int srcz,
					  float Tend)

{
    int i, j, nt, NR, NZ, NTMAX,  nzcc, nzww, nrsrc,js;
    double dx, dt, dum, muvisc;
    double *zcc, *ccv, *zww, *wwv, *zvec, *zz, *rvec, *rvecx, *rho;
    double **pn, **p, **vx, **vxn, **vy, **vyn;


    char infilename[80], ofilename[80], logfilename[80];
	time_t tstart, tend;
	FILE *flogid;
	
	sprintf(logfilename,"RunLogfile.txt"); flogid=fopen(logfilename,"w");
	
	tstart = time(0);
	fprintf(flogid, "Run started local time: %s \n", asctime(localtime(&tstart)));
	printf("Run started local time: %s \n", asctime(localtime(&tstart)));

// start timing
	timestamp();

//.........................................................................
// INPUT THE ENVIRONMENTAL 1D MODEL  - first TOPOGRAPHY

    double *rtopo, *topo;
    int nntop, itmin;
    sprintf(infilename,"topo.in");
    inputascii(infilename,rmax,rtopo,topo,nntop);
    double topomin = vector_min(topo,itmin,1,nntop);
    double topomax = vector_max(topo,itmin,1,nntop);
    printf("topomin + %lf,nntop = %d\n", topomin,nntop);
    if (rtopo[1] > 0 | rtopo[nntop] < rmax) {
        printf("FATAL ERROR: topography profile does not span 0-rmax");
        printf("\n rtopo[1]=%lf\n", rtopo[1]);
        printf(" rtopo[last]=%lf, rmax=%d\n",rtopo[nntop],rmax);
        fprintf(flogid,"FATAL ERROR: input topography does not span range");
        exit(-1);
    }
    if (topomin != zmin) {
        zmin = topomin;
        printf("WARNING: CODE HAS RESET ZMIN to = %d\n", zmin);
        fprintf(flogid, "RESET ZMIN to = %d due to topography \n", zmin);
    }
    if (topomax > zmax) {
        printf("FATAL ERROR: top of model must be > %lf m\n", topomax);
        exit(-1);
    }
    if (srcz<topo[1]) {
        srcz = topo[1];
        printf("WARNING: CODE HAS RESET SOURCE ALTITUDE to  %d\n", srcz);
        fprintf(flogid, "RESET SOURCE ALTITUDE to = %d due to topography \n", zmin);
    }
    
// input the 1D environmental model - sound speed
    
    sprintf(infilename,"sndprof.in");
    inputascii(infilename,zmax,zcc,ccv,nzcc);
	if (zcc[1] > zmin | zcc[nzcc] < zmax) {
		printf("FATAL ERROR: input sound profile does not span zmin-zmax");
		printf("\n zcc[1]=%lf, zmin=%d\n", zcc[1],zmin);
		printf(" zcc[last]=%lf, zmax=%d\n",zcc[nzcc],zmax);
		fprintf(flogid,"FATAL ERROR: input sound profile does not span zmin-zmax");
		exit(-1);
	}

// input the 1D environmental model - now wind speed
    
    sprintf(infilename,"wndprof.in");
    inputascii(infilename,zmax,zww,wwv,nzww);
	if (zww[1] > zmin | zww[nzww] < zmax) {
		printf("FATAL ERROR: input wind profile does not span zmin-zmax");
		printf("\n zww[1]=%lf, zmin=%d\n", zww[1],zmin);
		printf(" zww[last]=%lf, zmax=%d\n",zww[nzww],zmax);
		fprintf(flogid,"FATAL ERROR: input wind profile does not span zmin-zmax");
		exit(-1);
	}
    
//    muvisc = 0.000018*(4.0/3.0);    // =(4/3)*mu (no bulk viscosity)
	muvisc = 0.0;

//.........................................................................
    
// compute SPATIAL GRID DISCRETIZATION 
	int ndx = 10;
    dx = compute_dx(ccv,wwv,rmax,zmin,zmax,fsrc,NR,NZ,ndx);
    printf(" dx = %f, NR = %d, NZ = %d \n", dx,NR,NZ);
	fprintf(flogid, "\t dx = %f m,  NR = %d,  NZ = %d \n", dx,NR,NZ);
	fprintf(flogid, "\t viscosity = %f \n", muvisc);

//.........................................................................
// DISCRETIZE THE TOPOGRAPHY, WIND & SOUND SPEED PROFILES,
    
    rvec=zeros(NR); rvecx=zeros(NR); zvec=zeros(NZ); zz = zeros(NZ+1);
    model_grids(dx,NR,NZ,zmin,rvec,rvecx,zvec,zz);
    
    double *c1d=zeros(NZ), *w1d=zeros(NZ), *toporng = zeros(NR);
    toporng = linear_interp(rtopo,topo,rvec);
    sprintf(ofilename,"outputtopo.txt"); outputascii(ofilename,toporng);
    
    c1d = cubic_spline(zcc,ccv,zvec,0.,0.);
	sprintf(ofilename,"outputss.txt"); outputascii(ofilename,c1d);

//    printf("c1d \n");printvector(c1d,NZ); exit(1);
//    w1d = linear_interp(zww,wwv,zvec);
	w1d = cubic_spline(zww,wwv,zvec,0.,0.);
	sprintf(ofilename,"outputww.txt"); outputascii(ofilename,w1d);
 
//.........................................................................
// create ijrtopo vector, containing ij pairs where vx=0 (rigid boundary)
// & jztopo vector, containing j locations where vz=0 (rigid boundary)
    
    int **ijrtopo, nij;
    int *jztopo = izeros(NR), *jxtopo = izeros(NR+1);
    topoSetup(toporng,zz,jztopo,jxtopo);
    
    sprintf(ofilename,"outputjztopo.txt"); outputascii(ofilename,jztopo);
    sprintf(ofilename,"outputjxtopo.txt"); outputascii(ofilename,jxtopo);
    
    printf("Back from topoSetup, nij = %d\n",nij);
	
//.........................................................................
// CREATE RHO vector
	double *rhoprof;
	double rhomin = 0.00000003;
    rhoprof = compute_rhoprof(ccv,zcc,nzcc,grav); //layering due to gravity
	double rmin = vector_min(rhoprof,i,1,nzcc);
	
	if (rmin < rhomin) {
		printf(" rmin = %14.11f, rhomin = %14.11f \n",rmin,rhomin);
		fprintf(flogid,"density profile limited to rho >  %14.11f \n",rhomin);
		double *rhoint = fixrhoprof(zcc, rhoprof,rhomin);
//		sprintf(ofilename,"outputrhoi.txt"); outputascii(ofilename,rhoint);
		rho = cubic_spline(zcc,rhoint,zvec,0.,0.);
		delete[] rhoint;
	} else 	{
		rho = cubic_spline(zcc,rhoprof,zvec,0.,0.); 
	}
	sprintf(ofilename,"outputrho.txt"); outputascii(ofilename,rho);
    
    delete[] rhoprof; delete[] ccv; delete[] zcc;
    delete[] wwv; delete[] zww; delete[] zvec;
	
//.........................................................................
// compute TEMPORAL DISCRETIZATION given sound & wind speeds, and dx 
    
    dt = compute_dt(c1d,w1d,rho,Tend,fsrc,dx,muvisc,NTMAX);
    printf(" dt= %12.7f, Tfinal= %f s, NTMAX= %d\n", dt, NTMAX*dt, NTMAX);
	fprintf(flogid,"\t dt = %9.6f s, %d time steps to run to T = %f s\n", dt, NTMAX, NTMAX*dt);
	
//.........................................................................
// receiver locations 
// put receivers every rinc m starting with position right below source
    
    int rinc = 200, nrec = floor(rmax/rinc)-2;

	if (NZ>=100000) {
        printf(" FATAL ERROR \n");
        printf(" NZ = %d is >= 100000 \n",NZ);
		printf(" SUGGESTION: decrease fsrc_max or rmax & recompile \n");
		fprintf(flogid,"FATAL ERROR:  NZ %d is >= 100000 \n",NZ);
        exit(-1);  
    } 	
	
    int *irec = new int[nrec+1];       // receiver locs (radial indices)
	int *jrec = new int[nrec+1];       // receiver locs (altitude indices) 
    double **precvrs=zeros(nrec,NTMAX);     // P response at receivers
	printf("%d receivers, %d meters apart\n\n",nrec,rinc);
	fprintf(flogid,"\t %d receivers, %d meters apart\n",nrec,rinc);
    
    irec[0] = nrec;
    for (i=1; i<=nrec; i++) {			//print location @ end of precvrs
        irec[i] = 1+int(round(i*rinc/dx));
		jrec[i] = jztopo[irec[i]];
//        printf(" i= %d, irec = %d, at %lf   ",i,irec[i],irec[i]*dx-dx/2);
//        printf(" jrec = %d, at %lf \n",jrec[i],jrec[i]*dx-dx/2);
//        precvrs[i][0] = double(jrec[i]*100000+irec[i]);
        precvrs[i][0] = double(irec[i]) + double(jrec[i])/100000;
    }
	
//.........................................................................
// make multipliers for time the stepping loop (acoustic terms)
	
    double *Vxacterm = zeros(NZ); double *Vyacterm = zeros(NZ);
	double *Pacterm = zeros(NZ); double *pcylterm = zeros(NZ);
	double *Wterm = zeros(NZ); 	double *swterm = zeros(NZ);
	
    double **vxacoust = zeros(NR+1,NZ); double **vyacoust= zeros(NR,NZ+1);
    double **pacoust = zeros(NR,NZ);
	
	acoustSetup(Vxacterm,Vyacterm,Pacterm,pcylterm,Wterm,swterm,rho,c1d,
				w1d,dt,dx);
	
//.........................................................................
// create multipliers for the time stepping loop (viscosity terms)
	
    bool ifvisc = false;
	double *viscylr, *viscylz, *viscylrr, *viscx, *viscy, *muArt;
	double **vxvisc = zeros(NR+1,NZ);double **vyvisc= zeros(NR,NZ+1);

    if (muvisc!=0.0) {
		ifvisc = true;
	    viscx = zeros(NZ); viscy = zeros(NZ+1); muArt = zeros(NZ);
		viscylr = zeros(NZ); viscylrr = zeros(NR); viscylz = zeros(NZ+1);		
		printf(" computing viscosity terms \n");
// XXXXXX: does this need to be changed for wind???
		viscSetup(viscx,viscy,viscylr,viscylz,viscylrr,muArt,muvisc,rho,
					c1d,dt,dx,NR,NZ,fsrc);
		sprintf(ofilename,"muArt.txt"); outputascii(ofilename,muArt);
//		exit(1);
	}
	
//.........................................................................
// initialize Gaussian spatial source
    
    pn = zeros(NR,NZ); 
    nrsrc = newinit_source(NR,zmin,NZ,srcz,dx,pn,c1d,fsrc,jztopo,js);
    nt=0;
    sprintf(ofilename,"pr%2.2d.bin",0);
//    printf("for time step= %d, compute values to nrsrc = %d\n",0,nrsrc);
    outputBinary_2d(ofilename,pn,nt,dx,dt,zmin,nrsrc,js+nrsrc);
		
//.........................................................................
// initialize field variables and do the FIRST TIME STEP
    
    p = zeros(NR,NZ); vx = zeros(NR+1,NZ); vy = zeros(NR,NZ+1); 
	vxn = zeros(NR+1,NZ); vyn = zeros(NR,NZ+1);
	
	TstepOne(vx, vy, p, vxn, vyn, pn, Vxacterm, Vyacterm, Pacterm, 
			 pcylterm, rvec, jztopo, jxtopo, dt, dx, nrsrc);
//    printf(" calling applyTopo \n");
//    applyTopo(jztopo,vy, ijrtopo, vx, nij);
//    printf(" back from applyTopo \n");
	
//.........................................................................
// initialize the field variables and matrices for the boundary regions
    
    int iebc, jebc, iefbc, ibfbc, jbbc, ibbc;
    double **vxbcb, **vybcb, **pxbcb, **pybcb;  //  back PML bdy fields 
    double **vxbcr, **vybcr, **pxbcr, **pybcr;  //  right PML bdy fields
    
    double **cavybcb, **cbvybcb, **cavxbcb, **cbvxbcb; //back PML constants
    double **dapybcb, **dbpybcb, **dapxbcb, **dbpxbcb; //back PML constants
    double *cavy2=zeros(NR), *cbvy2=zeros(NR), *pcylbcb;  //back PML cons
    
    double **cavxbcr,**cbvxbcr,**dapxbcr,**dbpxbcr,**pcylbcr; //right PML
    double *cavx2=zeros(NZ), *cbvx2=zeros(NZ), *cwrit = zeros(NZ);
    
    iebc = 10; jebc=iebc;            //thickness of  PML regions
    iefbc=NR+iebc; ibfbc=iefbc+1; jbbc=jebc+1; ibbc=iebc+1;
    
    //fields in back and right boundaries
	double *cwtop = zeros(iefbc);
    vxbcb=zeros(ibfbc,jebc); vybcb=zeros(ibfbc,jbbc);
    pxbcb=zeros(iefbc,jebc); pybcb=zeros(iefbc,jebc);
	double **vxbcbn=zeros(ibfbc,jebc); double **vybcbn=zeros(ibfbc,jbbc);
    double **pxbcbn=zeros(iefbc,jebc); double **pybcbn=zeros(iefbc,jebc);
    vxbcr=zeros(ibbc,NZ); vybcr=zeros(iebc,NZ+1);
    pxbcr=zeros(iebc,NZ); pybcr=zeros(iebc,NZ);
	double **vxbcrn=zeros(ibbc,NZ); double **vybcrn=zeros(iebc,NZ+1);
    double **pxbcrn=zeros(iebc,NZ); double **pybcrn=zeros(iebc,NZ);
    
    //  constants in back PML boundary
    cavybcb=zeros(iefbc,jbbc); cbvybcb=zeros(iefbc,jbbc);
    dapybcb=zeros(iefbc,jebc); dbpybcb=zeros(iefbc,jebc);
    cavxbcb=zeros(iefbc,jebc); cbvxbcb=zeros(ibfbc,jebc);
    dapxbcb=zeros(iefbc,jebc); dbpxbcb=zeros(iefbc,jebc);
    pcylbcb=zeros(iefbc);
    
    //  constants in right PML boundary
    cavxbcr=zeros(ibbc,NZ); cbvxbcr=zeros(ibbc,NZ); dapxbcr=zeros(iebc,NZ);
    dbpxbcr=zeros(iebc,NZ); pcylbcr=zeros(iebc,NZ); 
	
// now allows for 1d wind
    PMLsetup(iebc,NR,NZ,dt,dx,c1d,w1d,rho,cavx2,cbvx2,cavy2,cbvy2,cwtop,cwrit,cavxbcb,
		cbvxbcb, cavybcb, cbvybcb, dapxbcb, dbpxbcb, dapybcb, dbpybcb,
		cavxbcr, cbvxbcr, dapxbcr, dbpxbcr,pcylbcb, pcylbcr);
    
// set up the vector for the chasing boundary condition, assume 1-D propagation
// uses 1st order accurate Mur 1D BC
// see ECE6340 Lecture 17.4 Mur boundary condition for 1D FDTD on youtube
    double *rbc = zeros(NZ);
    for (j=1; j<=NZ; j++) {					// main grid
        dum = c1d[j]-w1d[j];
        rbc[j] = (dum*dt-dx)/(dum*dt+dx);
    }
     
	double c1d1 = c1d[1];
	delete[] rho; delete[] c1d; delete[] w1d;

//.........................................................................

// time-stepping loop, spatial grid is staggered, not time grid
// Vx : vx[i=1][j]=0 by symmetry
// Vy : vy[i][j=1]=0 at the rigid surface

    timestamp();
	
	int ndump = 1000;
	double cleft = 295.0, cright = 356.0;
	int iwindoright = int(fmin(nrsrc+ndump*dt*cright/dx+100,NR));
	int jwindotop = int(fmin(js+nrsrc+ndump*dt*cright/dx+100,NZ));
	
	int iwindoleft = 2, i1=iwindoleft-1;
	double Tdelay = (zmax-2*zmin+srcz+nrsrc*dx)/cleft;
	printf("\t Tdelay = %lf, cleft=%lf \n", Tdelay, cleft);
	printf("\t iwindoright = %d, jwindotop = %d\n\n", iwindoright,jwindotop);

    fprintf(flogid,"\t Tdelay = %lf, cleft=%lf \n", Tdelay, cleft);
	fprintf(flogid,"\t iwindoright = %d, jwindotop = %d\n\n" ,iwindoright,jwindotop);
	int iwindRbdy = iwindoright;
	bool iright = false, ileft = true, jtop=false;
    
// put in a lightly absorbing column in front of chasing boundary
    int bufflen = 200;
    double habsorb, absorb = 0.006;
    habsorb = (1 - absorb/2);
    absorb = (1 - absorb);
	
	fprintf(flogid, "\n\n START of TIME STEP LOOP \n\n");
    
	tstart = time(0);
    for (nt=1; nt<=NTMAX; nt++) {
/*		swap2d(vxn,vx,1,iwindoright+1,1,jwindotop);
        swap2d(vyn,vy,1,iwindoright,1,jwindotop+1); 
        swap2d(pn,p,1,iwindoright,1,jwindotop);
 */
// use new code suggested by Jonathon Koenig
        std::swap(vxn,vx);
        std::swap(vyn,vy);
        std::swap(pn,p);
	
		acoustTerm(pacoust, vxacoust, vyacoust, pn, vxn, vyn, Vxacterm,
			Vyacterm, Pacterm, pcylterm, Wterm, swterm, rvec, iwindoleft, 
			iwindoright, jwindotop, ileft, iright);
        
        if (ileft) {
			
            for (i=i1; i<=iwindoright; i++) {
                for (j=1; j<=jwindotop; j++) {					// main grid
                    vx[i][j] = vx[i][j] + vxacoust[i][j];
                    vy[i][j] = vy[i][j] + vyacoust[i][j];
                    p[i][j] = p[i][j] + pacoust[i][j];
                }
            }
        
        } else {    // put in a slightly absorbing column
            
            for (i=i1; i<=i1+bufflen; i++) {
                for (j=1; j<=jwindotop; j++) {					
                    vx[i][j] = vx[i][j] + vxacoust[i][j];
                    vy[i][j] = vy[i][j] + vyacoust[i][j];
                    p[i][j] = absorb*p[i][j] + habsorb*pacoust[i][j];
                }
            }
            
            i = i1; // overwright 1st column of vx with 1st order Mur 1D BC
            for (j=1; j<=jwindotop; j++) {
                vx[i][j] = vxn[i+1][j]+rbc[j]*(vx[i+1][j]-vxn[i][j]);
            }
            
            for (i=i1+bufflen+1; i<=iwindoright; i++) {
                for (j=1; j<=jwindotop; j++) {					// main grid
                    vx[i][j] = vx[i][j] + vxacoust[i][j];
                    vy[i][j] = vy[i][j] + vyacoust[i][j];
                    p[i][j] = p[i][j] + pacoust[i][j];
                }
            }
        }
        
        applyTopo(jztopo,vy, jxtopo, vx);     //set V=0 within rigid boundary


        if (ifvisc) {	
			viscxterm(vxvisc,  vx, viscx, viscylr, viscylrr, rvecx, 
					  iwindoleft, iwindoright,jwindotop);
			visczterm(vyvisc, vy, viscy, viscylz, rvec,  iwindoleft, 
					  iwindoright,jwindotop);	 
            
            for (i=iwindoleft; i<=iwindoright; i++) {				
                for (j=1; j<=jwindotop; j++) {					// main grid
					vx[i][j] = vx[i][j] + vxvisc[i][j];
					vy[i][j] = vy[i][j] + vyvisc[i][j];
                }
            }
            
// do I need to do absorbing colum here too (or here instead of above?)
			
			if (ileft){						// source region vx[1][j]=0
                for (j=1; j<=jwindotop; j++) vy[1][j] = vy[1][j] + vyvisc[1][j];
            }
            
            applyTopo(jztopo,vy, jxtopo, vx);   //set V=0 within rigid boundary
			
        } 
		
// save pressure response at receivers
        for (j=1; j<=nrec; j++) precvrs[j][nt] = p[irec[j]][jrec[j]];

//.........................................................................
// do PML updates for boundaries
		if (iright | jtop) {
			topBound(iwindoright, iwindRbdy, jebc, NZ, vxbcb, vxbcbn, vybcb, 
				vybcbn, pxbcb, pxbcbn, pybcb, pybcbn, vy, pn, vyn, cavy2,
				cbvy2, cwtop, cavxbcb,  cbvxbcb, cavybcb, cbvybcb, dapxbcb, 
				dbpxbcb, dapybcb, dbpybcb, pcylbcb);		
			
// must call topBound before rightBound (because of swaps)
			if (iright) {	// only call rightBound after Tdelay seconds	
//.........................................................................
				rightBound(iebc, NZ, NR, Vyacterm, vxbcr, vxbcrn, vybcr, vybcrn,
					pxbcr, pxbcrn, pybcr, pybcrn, pybcb, pybcbn, pxbcbn,vybcbn,
					vx, vxn, pn, cavx2, cbvx2, cwrit,cavxbcr, cbvxbcr, dapybcb, 
					dbpybcb, dapxbcr, dbpxbcr, Pacterm, pcylterm, rvec);	
			}
		}
// end of boundary condition updates
//.........................................................................
        
        if ( nt % ndump ==0 ) {           //print interim results
// output pressure		
            sprintf(ofilename,"pr%2.2d.bin",nt);
			printf("nt = %d, filename %s \n",nt,ofilename);
			outputBinary_2d(ofilename,p,nt,dx,dt,zmin,iwindoright,jwindotop); 
			
// output velocities (vxn + next time step, vx)
/*
            sprintf(ofilename,"vx%2.2d.bin",nt);
            outputBinary_2d(ofilename,vx,nt,dx,dt,zmin,iwindoright+1,jwindotop);
            sprintf(ofilename,"vy%2.2d.bin",nt);
			outputBinary_2d(ofilename,vy,nt,dx,dt,zmin,iwindoright,jwindotop+1); 
			
            sprintf(ofilename,"vxn%2.2d.bbin",nt);
            outputBinary_2d(ofilename,vxn,nt,dx,dt,zmin,iwindoright+1,jwindotop);
            sprintf(ofilename,"vyn%2.2d.bbin",nt);
			outputBinary_2d(ofilename,vyn,nt,dx,dt,zmin,iwindoright,jwindotop+1);
			
			if (ifvisc) {
// output viscosity term (in velocity) (either stabilized or not)
            sprintf(ofilename,"vxmuvisc%2.2d.bbin",nt);
            outputBinary_2d(ofilename,vxvisc,nt,dx,dt,zmin,iwindoright,jwindotop);
            sprintf(ofilename,"vymuvisc%2.2d.bbin",nt);
			outputBinary_2d(ofilename,vyvisc,nt,dx,dt,zmin,iwindoright,jwindotop);
			}
			
// output acoustic term (in velocity)
            sprintf(ofilename,"vxacoust%2.2d.bbin",nt);
            outputBinary_2d(ofilename,vxacoust,nt,dx,dt,zmin,iwindoright,jwindotop);
            sprintf(ofilename,"vyacoust%2.2d.bbin",nt);
			outputBinary_2d(ofilename,vyacoust,nt,dx,dt,zmin,iwindoright,jwindotop); 
*/ 					
// output interim response at receivers 
			int nlive = int(fmin(nrec,ceil(nt*dt*c1d1/rinc)));
			printf(" printing out response at %d receivers\n",nlive);
            sprintf(ofilename,"prespCyl%2.2d.bin",nt);
            outputBinary_2d(ofilename,precvrs,NTMAX,dx,dt,zmin,nlive,nt);

			iwindoright = int(fmin(nrsrc+(nt+ndump)*dt*cright/dx+100,NR));
			jwindotop = int(fmin(js+nrsrc+(nt+ndump)*dt*cright/dx+100,NZ));
			iwindRbdy = iwindoright;			
			if (iwindoright == NR)	{
				iright=true;
				iwindRbdy = iefbc;
			}   

/* XXXX more generally, should set the left window a given distance
        behind the right window based on:
 a) how long the maximum coda duration is expected to be. ie for 30s coda 
        - make window about 10 km wide + buffers at each end
 b) the height of the window. The sound field is nearly spherical near 
        the source so it should be about as wide as it is high, at least
        near the source (sound field steepens further from source)
 (idea: could "chase" at same velocity as right window)
 */
			iwindoleft = int(fmax(2, (nt*dt-Tdelay)*cleft/dx));
			i1=iwindoleft-1;
			if (iwindoleft>= 5)	{
                cleft = 325;
				ileft=false;
				printf("zero fields at left \n");
// boundaries too? (vxbcb, vxbcbn, vybcb, vybcbn, pxbcb, pxbcbn, pybcb, pybcbn)
                for (i=1; i<=i1-2; i++) {
                    for (j=1; j<=NZ; j++) {					// main grid
                        vx[i][j] = 0;  vy[i][j] = 0;  p[i][j] = 0;
                        vxn[i][j] = 0; vyn[i][j] = 0; pn[i][j] = 0;
                    }
                    vy[i][NZ+1] = 0; vyn[i][NZ+1] = 0;
                }
            }
            
			if (jwindotop == NZ) {
				jtop = true;
			}
			printf("left & right edges @ i=%d, %d\n",iwindoleft,iwindoright);
			printf("right bdy index=%d\n",iwindRbdy);
			printf("top edge @ j= %d\n",jwindotop);
									
			fprintf(flogid,"\n\n Interim output at nt=%d steps, %f s\n", nt,nt*dt);
			
			tend = time(0);
			fprintf(flogid, " Local time : %s ", asctime(localtime(&tend)));
			fprintf(flogid, " Time for %d steps: %g s\n",nt,difftime(tend,tstart));	
			fprintf(flogid, "\t iwindoleft = %d, iwindoright = %d\n",iwindoleft,iwindoright);
			fprintf(flogid, "\t jwindotop = %d\n",jwindotop);
			printf("Time for %d steps: %g s\n",nt,difftime(tend,tstart));			
			timestamp();
        }				// end printing interim solution		
    }                                   // end time-stepping loop
    timestamp();
    
    // print final responses
    sprintf(ofilename,"prespCyl%2.2d.bin",NTMAX);
    outputBinary_2d(ofilename,precvrs,NTMAX,dx,dt,zmin);
	tend = time(0);
	fprintf(flogid, "\n\n Finished at local time : %s", asctime(localtime(&tend)));
	fclose(flogid);

    return;
}
//#########################################################################