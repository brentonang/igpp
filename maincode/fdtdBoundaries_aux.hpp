//#########################################################################

// functions having to do with the boundary layers at top and right

//#########################################################################

void PMLsetup(int iebc, int NR,int NZ,double dt,double dx,double *cc,
		double *ww, double *rho,double *cavx2,double *cbvx2,double *cavy2,
		double *cbvy2, double *cwtop, double *cwrit, double **cavxbcb,double **cbvxbcb,double **cavybcb,
		double **cbvybcb,double **dapxbcb,double **dbpxbcb,double **dapybcb,
		double **dbpybcb,double **cavxbcr,double **cbvxbcr,double **dapxbcr,
		double **dbpxbcr,double *pcylbcb, double **pcylbcr);

//######################################################################### 

void topBound(int iwindoright, int iwindRbdy, int jebc, int NZ,  
         double **vxbcb, double **vxbcbn, double **vybcb, double **vybcbn,
         double **pxbcb, double **pxbcbn, double **pybcb, double **pybcbn,
         double **vy, double **pn,double **vyn,double *cavy2,double *cbvy2,
         double *cwtop,double **cavxbcb, double **cbvxbcb,double **cavybcb, 
         double **cbvybcb, double **dapxbcb, double **dbpxbcb, 
         double **dapybcb, double **dbpybcb, double *pcylbcb);

//#########################################################################  

void rightBound(int iebc, int NZ, int NR, double *dtdxorz,
		double **vxbcr, double **vxbcrn, double **vybcr, double **vybcrn,
		double **pxbcr, double **pxbcrn, double **pybcr, double **pybcrn,
		double **pybcb, double **pybcbn, double **pxbcbn, double **vybcbn,
		double **vx, double **vxn, double **pn,double *cavx2,double *cbvx2,
		double *cwrit, double **cavxbcr, double **cbvxbcr, double **dapybcb, 
		double **dbpybcb, double **dapxbcr, double **dbpxbcr, 
		double *dtdxrhocsq, double *pcylterm, double *rvec);

//#########################################################################

void topoSetup(double toporng[], double zvec[], int jztopo[],int jxtopo[]);

//int **topoSetup(double toporng[], double zvec[], int jztopo[], int &nij);
/*  PURPOSE:  set up indices to enforce the rigid boundary conditions at
    a rigid boundary. Particle velocities are zero at a rigid boundary.
    the size of the jztopo vector, which contains j locations where vz=0 
    is known in advance. The ijrtopo vector gives i,j paris where vertical
    boundaries are zero. Its size is not known in advance, and has to be
    set within the program (it is small for a nearly flat boundary, large
    for a very bumpy surface. 
    the vz and vx values at the boundary are set to zero in applyTopo
    INPUT: 
        toporng - topography as a function of range
        zvec - altitude vector
    OUTPUT:
        jztopo - vector of j indices where vz is set to zero in applyTopo
 NO       int ** output: 2 X NIJ array of i,j indices where vx is set to 0
 YES    jxtopo - vector of j indices, vx is !0 @ all points >=jxtopo
*/

//#########################################################################

void applyTopo(int *jztopo, double **vz, int *jxtopo, double **vx);
// void applyTopo(int *jztopo, double **vz, int **ijpairs, double **vx,
//                int nij);
/*  PURPOSE:  uses indices computed in topoSetup to enforce rigid boundary
            condition at topography
 INPUT:
    jztopo - vector of j indices where vz is set to zero in applyTopo
NO    ijpairs -  2 X NIJ array of i,j indices where vx is set to 0
    jxtopo - vector of j indices  vx !=0  for j>=jxtopo
NO    nij - length of ijpairs
    vx, vz - particle velocity arrays, altered by this function
*/

//#########################################################################