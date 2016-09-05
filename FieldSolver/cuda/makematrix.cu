/*

makematrix.cu

Constructs matrix equation A x = b for local minimization problem 
Exports matrix A and vector b into plain text

Written by Hee Sok Chung at ANL
July 10, 2016

Modified by Ran Hong for cuda compatibility

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>


/* Legendre functions are defined without normalization factors in order to 
 * avoid overflow from gamma function */

/* hypergeometric function 1F2 */
double h1f2(double a, double b, double c, double z);

/* LegendreQ(m-1/2,n,cosh(x)) */
double LegendreQ(double m, double n, double z);

/* derivative of LegendreQ(m-1/2,n,cosh(x)) */
double DLegendreQ(int m, int n, double z);

/* coordinate transformation r, z -> zeta */
double zetaf(double rho, double z, double r0);

/* coordinate transformation r, z -> eta */
double etaf(double rho, double z, double r0);

/* main subroutine that calculates matrix and vector elements */
void cpx(double **datain, 
    double *cc, double *cs, double *sc, double *ss,
    double **matrix, 
    double *zeta, double *eta, 
    int **ccidx, int ndata, int dim1, 
    double zt0, double rr, int *sliceflag);

/**
 * CUDA Kernel Device code
 *
 * Computes the vector addition of A and B into C. The 3 vectors have the same
 * number of elements numElements.
 */
__global__ void vectorProd(const double *d_tcc, const double *d_tcs, const double *d_tsc, const double *d_tss, double *d_M, int dim_v)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;
  int ii=0;
  int jj=0;
  int dim_M = 4*dim_v+2;
  double local_M[4][4];


  if (i < dim_v && j< dim_v){
    for (ii=0;ii<4;ii++){
      for (jj=0;jj<4;jj++){
	local_M[ii][jj]=d_M[(4*i+ii+1)*dim_M+4*j+jj+1];
      }
    }

    local_M[0][0]+=d_tcc[i]*d_tcc[j];
    local_M[1][0]+=d_tcs[i]*d_tcc[j];
    local_M[2][0]+=d_tsc[i]*d_tcc[j];
    local_M[3][0]+=d_tss[i]*d_tcc[j]; 

    local_M[0][1]+=d_tcc[i]*d_tcs[j];
    local_M[1][1]+=d_tcs[i]*d_tcs[j];
    local_M[2][1]+=d_tsc[i]*d_tcs[j]; 
    local_M[3][1]+=d_tss[i]*d_tcs[j];

    local_M[0][2]+=d_tcc[i]*d_tsc[j]; 
    local_M[1][2]+=d_tcs[i]*d_tsc[j]; 
    local_M[2][2]+=d_tsc[i]*d_tsc[j]; 
    local_M[3][2]+=d_tss[i]*d_tsc[j]; 

    local_M[0][3]+=d_tcc[i]*d_tss[j]; 
    local_M[1][3]+=d_tcs[i]*d_tss[j]; 
    local_M[2][3]+=d_tsc[i]*d_tss[j]; 
    local_M[3][3]+=d_tss[i]*d_tss[j]; 

    for (ii=0;ii<4;ii++){
      for (jj=0;jj<4;jj++){
	d_M[(4*i+ii+1)*dim_M+4*j+jj+1]=local_M[ii][jj];
      }
    }
  }
}


int main (int argc, char **argv) {
  int nmax, mmax, ndata, dim1, dim2x, vecdim;
  //  Max. number of harmonics in azimuthal direction 
  nmax=200; 
  //  Max. number of harmonics in poloidal direction 
  mmax=8; 
  //  Matrix and vector dimensions 
  dim1=(nmax+1)*(mmax+1);
  dim2x=4*dim1;
  double rr;
  double *cc, *cs, *sc, *ss;
  double **matrix, **datain;
  double datatmp;
  double zeta[30], eta[30];
  double zt0, bzero;
  double angles[30], radii[30];
  int **ccidx;
  int **mnidx;
  int *sinz, *vecidx;
  int i, j, id, md, nd, prb ;
  int vecc;
  FILE *input;
  FILE *inputc;
  //    procnum=11; // number of thread to use 
  //    printf("Using %d threads\n",procnum);

  //  declare vector and matrix elements as pointers 
  cc=(double *)malloc((dim1+2)*sizeof(double));
  cs=(double *)malloc((dim1+2)*sizeof(double));
  sc=(double *)malloc((dim1+2)*sizeof(double));
  ss=(double *)malloc((dim1+2)*sizeof(double));

  matrix=(double **)malloc((dim2x+2)*sizeof(double*));
  for (i=0;i<(dim2x+2);i++) {
    matrix[i]=(double *)malloc((dim2x+2)*sizeof(double));
  }
  sinz=(int *)malloc((dim2x+2)*sizeof(int)); // for identifying exact zeros 
  vecidx=(int *)malloc((dim2x+2)*sizeof(int)); // vector index 
  ccidx=(int **)malloc((dim1+2)*sizeof(int*)); // matrix single index -> double index 
  for (i=0;i<(dim1+2);i++) {
    ccidx[i]=(int *)malloc(3*sizeof(int));
  }
  mnidx=(int **)malloc((dim1+2)*sizeof(int*)); // matrix double index -> single index
  for (i=0;i<(dim1+2);i++) {
    mnidx[i]=(int *)malloc((dim1+2)*sizeof(int));
  }

  //  count number of azimuthal slices 
  ndata=0;
  //  inputc=fopen("data40.txt", "r"); // open data for counting
  inputc=fopen("data52.txt", "r"); // open data for counting
  while(fscanf(inputc, "%lf", &datatmp)>0) {
    ndata++;
  }
  ndata=ndata/26; // azimuthal angle + 25 probes 
  fclose(inputc); // close data file used for counting 
  int sliceflag[ndata+1];
  printf("number of data = %d\n", ndata);
  //  this pointer will be used for data taking 
  datain=(double **)malloc((ndata+1)*sizeof(double*));
  for (id=0;id<=ndata;id++) {
    datain[id]=(double *)malloc((25+2)*sizeof(double));
  }
  //  distribute azimuthal slices to threads
  //    printf("%d\n",procnum);
  //    printf("%d\n",ndatap);
  i=1;
  for (id=1;id<=ndata;id++){
    sliceflag[id]=1;
  }
  //    printf("Max i = %d\n", i);

  bzero=61.789; // reference B field, chosen close to average B field 
  zt0=6.5; // toroidal harmonics will be normalized at zeta=zt0
  //  construct probe positions 
  angles[1]=0;
  radii[1]=0;
  for (i=2;i<=9;i++) {
    angles[i]=(i-2.)*M_PI/4;
    radii[i]=22.5;
  }
  for (i=10;i<=25;i++) {
    angles[i]=(i-10.)*M_PI/8;
    radii[i]=45.;
  }
  /*    for (i=1;i<=25;i++) {
	printf("probe %d location r=%lf, theta=%lf \n",i+1,radii[i],angles[i]);
	} */

  //  calculate toroidal coordinates 
  rr=7111.5; // toroid center. MUST NEVER COINCIDE WITH ACTUAL PROBE POSITION 
  printf ("\nComputing coordinates\n");
  for (i=1;i<=25;i++) {
    zeta[i]=zetaf(radii[i]*sin(angles[i])+7112., \
	-radii[i]*cos(angles[i]),rr);
    eta[i]=etaf(radii[i]*sin(angles[i])+7112., \
	-radii[i]*cos(angles[i]),rr);
  }
  /* initialize arrays */
  i=0;
  nd=0;
  md=0;
  printf ("\nInitializing coefficients\n");
  for (i=0;i<=dim2x;i++) {
    sinz[i]=1;
  }
  i=0;
  for (nd=0;nd<=nmax;nd++) {
    for (md=0;md<=mmax;md++) {
      i+=1;
      ccidx[i][1]=md; ccidx[i][2]=nd; 
      mnidx[md][nd]=4*(i-1);
      cc[i]=0.; 
      cs[i]=0.; 
      sc[i]=0.; 
      ss[i]=0.;
      if(nd==0) { // identify sin(0) = 0 due to n = 0
	sinz[4*(i-1)+2]=0;
	sinz[4*(i-1)+4]=0;
      }
      if(md==0) { // identify sin(0) = 0 due to m = 0
	sinz[4*(i-1)+3]=0;
	sinz[4*(i-1)+4]=0;
      }
    }
  }
  printf ("\nInitializing matrix\n");
  for (i=0;i<=dim2x;i++) {
    for (j=0;j<=dim2x;j++) {
      matrix[i][j]=0.;
    }
  }
  //    input =fopen("data40.txt", "r"); // this one's for actually reading data
  input =fopen("data52.txt", "r"); // this one's for actually reading data
  for (id=1;id<=ndata;id++) {
    fscanf (input, "%lf", &datatmp); // read azimuthal angle in degrees 
    datain[id][0]=datatmp;
    for (prb=1;prb<=25;prb++){ // loop over 25 probes 
      fscanf (input, "%lf", &datatmp); // read B-fields in kHz 
      datatmp=datatmp*0.001+61.7400000-bzero; // convert to MHz and offset
      datain[id][prb]=datatmp;
    }
  }
  fclose(input); // close data 

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  printf ("\nComputing matrix elements\n");
  //  Calculate matrix and vector elements 
  cpx(datain, cc, cs, sc, ss, matrix, 
      zeta, eta, ccidx, ndata, dim1, zt0, rr, sliceflag);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  for(i=1;i<=ndata;i++) {
    if (sliceflag[i]!=0) printf("slice %d status = %d\n", i, sliceflag[i]);
  }
  //  combine results from individual threads 

  FILE *out0, *out1, *outn, *outd;
  out0 = fopen ("outmatrix.dat", "w"); // matrix elements 
  out1 = fopen ("outvector.dat", "w"); // vector elements 
  outn = fopen ("outnames.txt", "w");  // index dictionary 
  outd = fopen ("vecdim.txt", "w");    // matrix dimensions

  fprintf(outn, "%d\n", mmax);
  fprintf(outn, "%d\n", nmax);

  for (i=1;i<=dim1;i++) {
    fprintf(outn, "1\t%d\t%d\n",ccidx[i][1], ccidx[i][2]);
    if(sinz[4*(i-1)+2]!=0) {
      fprintf(outn, "2\t%d\t%d\n",ccidx[i][1], ccidx[i][2]);
    }
    if(sinz[4*(i-1)+3]!=0) {
      fprintf(outn, "3\t%d\t%d\n",ccidx[i][1], ccidx[i][2]);
    }
    if(sinz[4*(i-1)+4]!=0) {
      fprintf(outn, "4\t%d\t%d\n",ccidx[i][1], ccidx[i][2]);
    }
  }
  vecdim=0;
  for (i=1;i<=dim1;i++) {
    vecdim++;
    //    fprintf(out1, "%.17g\n", cc[i]);
    fwrite(&(cc[i]),sizeof(double),1,out1);
    if(sinz[4*(i-1)+2]!=0) {
      vecdim++;
      //      fprintf(out1, "%.17g\n", cs[i]);
      fwrite(&(cs[i]),sizeof(double),1,out1);
    }
    if(sinz[4*(i-1)+3]!=0) {
      vecdim++;
      //      fprintf(out1, "%.17g\n", sc[i]);
      fwrite(&(sc[i]),sizeof(double),1,out1);
    }
    if(sinz[4*(i-1)+4]!=0) {
      vecdim++;
      //      fprintf(out1, "%.17g\n", ss[i]);
      fwrite(&(ss[i]),sizeof(double),1,out1);
    }
  }
  id=0;
  for (i=1;i<=dim2x;i++) {
    if(sinz[i]!=0) {
      id++;
      vecidx[id]=i;
    }
  }
  vecc=0;
  for (i=1;i<=vecdim;i++) {
    for (j=1;j<=vecdim;j++) {
      vecc++;
//      fprintf(out0, "%.17g\n", matrix[vecidx[i]][vecidx[j]]);
      fwrite(&(matrix[vecidx[i]][vecidx[j]]),sizeof(double),1,out0);
    }
  }
  fprintf(outd, "%d\n", vecdim);

  fclose(out0);
  fclose(out1);
  fclose(outn);
  fclose(outd);
  free(cc);
  free(cs);
  free(sc);
  free(ss);
  free(ccidx);
  return 0;
}

/* hypergeometric function 1F2 
   evaluated by truncating an infinite sum */
double h1f2(double a, double b, double c, double z){
  int i, imax;
  double err, tol, si, sii, f1f2;
  double errabs, errrel;
  imax=100000000; // 10^8 maximum iterations
  sii=1.; // initial term 
  err=1.; // estimated uncertainty 
  f1f2=1.; // initial contribution 
  i=0;  // iterator
  tol=pow(10.,-16.); // error tolerance 
  if (z<tol||z>1.-tol) { // z out of range or dangerously close to 0 or 1
    i=imax+1;
    f1f2=1./i; // to return inf or nan
    i=imax+1; // stop evaluation
  }
  while (err>tol&&i<=imax) {
    i++;
    si=(a+i-1.)*(b+i-1.)/(c+i-1.)*z/i*sii; // next term 
    f1f2+=si; //next term added 
    errabs=fabsl(si*z/(1.-z)); // estimated absolute uncertainty 
    errrel=fabsl(errabs/f1f2); // estimated relative uncertainty 
    if (errabs>errrel) { // choose larger one as error
      err=errabs;
    }
    else {
      err=errrel;
    }
    sii=si;
  }
  return f1f2;
}

/* Legendre function of the second kind, Q_{m+1/2}^n (cosh(z)).
   Normalized to remove gamma function.
   Regular inside the torus */
double LegendreQ(double m, double n, double z){
  double lq;
  lq=pow(tanh(z),n)/pow(cosh(z),m+.5);
  lq=lq*h1f2(.5*(m+n+.5),.5*(m+n+1.5),m+1.,1./cosh(z)/cosh(z));
  return lq;
}

/* Derivative of Q_{m+1/2}^n (cosh(z)).
   Normalized to remove gamma function.
 */
double DLegendreQ(int m, int n, double z){
  double dlq, lq1, lq2;
  if (m==0) {
    dlq=-1/(8*pow(cosh(z),1.5))/sinh(z);
    dlq=dlq*pow(tanh(z),n);
    dlq=dlq*( (4.*pow(sinh(z),2.)-8.*n)*                        \
	h1f2(n/2.+.25,n/2.+.75,1.,1/cosh(z)/cosh(z)) +   \
	(4.*n*n+8.*n+3.)*tanh(z)*tanh(z)*                \
	h1f2(n/2.+1.25,n/2.+1.75,2.,1/cosh(z)/cosh(z)));
  }
  else {
    lq1=pow(tanh(z),n)/pow(cosh(z),m+.5);
    lq1=lq1*h1f2(.5*(m+n+.5),.5*(m+n+1.5),m+1.,1./cosh(z)/cosh(z));
    lq2=pow(tanh(z),n)/pow(cosh(z),m-.5);
    lq2=lq2*h1f2(.5*(m+n-.5),.5*(m+n+.5),m*1.,1./cosh(z)/cosh(z));
    dlq=(m-.5)/tanh(z)*lq1-(2.*m)/sinh(z)*lq2;
  }
  return dlq;
}

/* coordinate transformation r, z -> zeta */
double zetaf(double rho, double z, double r0){
  double zetax;
  zetax=atanh(2.*rho*r0/(rho*rho+r0*r0+z*z));
  return zetax;
}

/* coordinate transformation r, z -> eta */
double etaf(double rho, double z, double r0){
  double etax, xx;
  int i;
  xx=2.*r0*z/(rho*rho-r0*r0+z*z);
  if (fabsl(xx)<0.001) {
    etax=1.;
    for (i=1;i<=10;i++) {
      etax=etax+pow(xx,2.*i)/(2.*i+1.)*cos(M_PI*i);
    }
    etax=etax*xx;
  }
  else {
    etax=atanl(xx);
  }
  if (rho<sqrtl(r0*r0-z*z)) etax=etax+M_PI;
  return etax;
}

void cpx(double **datain, 
    double *cc, double *cs, double *sc, double *ss,
    double **matrix, 
    double *zeta, double *eta, 
    int **ccidx, int ndata, int dim1, 
    double zt0, double rr, int *sliceflag) {
  double datatmp, zt, et, wgt, phi;
  double *tcc, *tcs, *tsc, *tss;
  tcc=(double *)malloc(sizeof*tcc*(dim1));
  tcs=(double *)malloc(sizeof*tcs*(dim1));
  tsc=(double *)malloc(sizeof*tsc*(dim1));
  tss=(double *)malloc(sizeof*tss*(dim1));
  int md, nd, i, id, prb, j ;
  double legq, dleq;
  cudaError_t err = cudaSuccess;

  //Calculate legq and dleq table
  double **Tlegq;
  double **Tdleq;
  Tlegq=(double **)malloc((dim1+1)*sizeof(double*));
  Tdleq=(double **)malloc((dim1+1)*sizeof(double*));
  for (int il=0;il<=dim1;il++) {
    Tlegq[il]=(double *)malloc(26*sizeof(double));
    Tdleq[il]=(double *)malloc(26*sizeof(double));
  }
  for (prb=1;prb<=25;prb++){ // loop over 25 probes 
    zt=zeta[prb];  // zeta coordinate at probe 
    for (i=1;i<=dim1;i++){
      md=ccidx[i][1]; // m index 
      nd=ccidx[i][2]; // n index 
      /* LegendreQ at probe */
      Tlegq[i][prb]=LegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
      /* Derivative of LegendreQ at probe */
      Tdleq[i][prb]=DLegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
    }
  }
  //Table calculation complete

  //Re-allocate matrix
  double *h_M = NULL;
  int dimM = 4*dim1+2;
  size_t sizeM = dimM*dimM*sizeof(double);
  h_M = (double *)malloc(sizeM);

  //device Malloc
  double *d_tcc = NULL;
  double *d_tcs = NULL;
  double *d_tsc = NULL;
  double *d_tss = NULL;
  int d_size = sizeof(double)*dim1;
  double *d_M = NULL;
  err = cudaMalloc((void **)&d_tcc, d_size);
  err = cudaMalloc((void **)&d_tcs, d_size);
  err = cudaMalloc((void **)&d_tsc, d_size);
  err = cudaMalloc((void **)&d_tss, d_size);
  err = cudaMalloc((void **)&d_M, sizeM);

  //copy matrix memory to device
  err = cudaMemcpy(d_M, h_M, sizeM, cudaMemcpyHostToDevice);

  for (id=1;id<=ndata;id++) { // loop over azimuthal slices 
    if (id!=0) { 
      sliceflag[id]=0;
      datatmp=datain[id][0];
      phi=datatmp*M_PI/180.;
      for (prb=1;prb<=25;prb++){ // loop over 25 probes 
	datatmp=datain[id][prb];
	zt=zeta[prb];  // zeta coordinate at probe 
	et=eta[prb];  // eta coordinate at probe 
	wgt=sqrt(cosh(zt)-cos(et)); // weight func in toroidal coordinates 
	//
	double Sinh_zt = sinh(zt);
	double Sin_et = sin(et);
	double Cosh_zt = cosh(zt);
	double Cos_et = cos(et);
	//
	for (i=1;i<=dim1;i++){
	  md=ccidx[i][1]; // m index 
	  nd=ccidx[i][2]; // n index 
	  /* LegendreQ at probe */
	  //legq=LegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
	  legq=Tlegq[i][prb];
	  /* Derivative of LegendreQ at probe */
	  //dleq=DLegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
	  dleq=Tdleq[i][prb];
	  // CC(m,n) coefficient 
	  tcc[i-1]=-Sinh_zt*Sin_et/rr*cos(nd*phi)*cos(md*et)*  \
		   (Sinh_zt/2./wgt*legq+wgt*dleq) \
		   -(1.-Cosh_zt*Cos_et)/rr* \
		   (-md*cos(nd*phi)*sin(md*et)*wgt*legq \
		    +cos(nd*phi)*cos(md*et)*Sin_et/2./wgt*legq);
	  // CS(m,n) coefficient 
	  tcs[i-1]=-Sinh_zt*Sin_et/rr*sin(nd*phi)*cos(md*et)* \
		   (Sinh_zt/2./wgt*legq+wgt*dleq) \
		   -(1.-Cosh_zt*Cos_et)/rr* \
		   (-md*sin(nd*phi)*sin(md*et)*wgt*legq \
		    +sin(nd*phi)*cos(md*et)*Sin_et/2./wgt*legq);
	  // SC(m,n) coefficient 
	  tsc[i-1]=-Sinh_zt*Sin_et/rr*cos(nd*phi)*sin(md*et)* \
		   (Sinh_zt/2./wgt*legq+wgt*dleq) \
		   -(1.-Cosh_zt*Cos_et)/rr* \
		   (md*cos(nd*phi)*cos(md*et)*wgt*legq \
		    +cos(nd*phi)*sin(md*et)*Sin_et/2./wgt*legq);
	  // SS(m,n) coefficient 
	  tss[i-1]=-Sinh_zt*Sin_et/rr*sin(nd*phi)*sin(md*et)* \
		   (Sinh_zt/2./wgt*legq+wgt*dleq) \
		   -(1.-Cosh_zt*Cos_et)/rr* \
		   (md*sin(nd*phi)*cos(md*et)*wgt*legq \
		    +sin(nd*phi)*sin(md*et)*Sin_et/2./wgt*legq);
	} // define tcc loop end

	//Start parallel computation on gpu
	//copy memory
	err = cudaMemcpy(d_tcc, tcc, d_size, cudaMemcpyHostToDevice);
	err = cudaMemcpy(d_tcs, tcs, d_size, cudaMemcpyHostToDevice);
	err = cudaMemcpy(d_tsc, tsc, d_size, cudaMemcpyHostToDevice);
	err = cudaMemcpy(d_tss, tss, d_size, cudaMemcpyHostToDevice);

	dim3 DimBlock (16,16);
	dim3 DimGrid (dim1/16+1, dim1/16+1);
//	printf("CUDA kernel launch with %d blocks of %d threads\n", DimGrid.x * DimGrid.y, 256);

	vectorProd<<<DimGrid, DimBlock>>>(d_tcc, d_tcs, d_tsc, d_tss, d_M, dim1);
	err = cudaGetLastError();

	if (err != cudaSuccess)
	{
	  fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
	  exit(EXIT_FAILURE);
	}

	for (i=1;i<=dim1;i++){
	  cc[i]+=tcc[i-1]*datatmp;
	  cs[i]+=tcs[i-1]*datatmp;
	  sc[i]+=tsc[i-1]*datatmp;
	  ss[i]+=tss[i-1]*datatmp;
	} // define vector loop end 
	// for each data end
      } // probe loop end 
      //    printf ("\nSlice %d (phi=%lf) done\n", id, phi);
      printf ("\nSlice %d of %d done\n", id, ndata);
    }
  } // id end 

  //copy matrix from device to host
  err = cudaMemcpy(h_M, d_M, sizeM, cudaMemcpyDeviceToHost);

  //copy matrix back to main
  for (i=0;i<dimM;i++){
    for (j=0;j<dimM;j++) {
      matrix[i][j]=h_M[i*dimM+j];
//      printf("%d,%d,%.17g\n",i,j,matrix[i][j]);
    }
  }
  // Free device global memory
  err = cudaFree(d_tcc);
  err = cudaFree(d_tcs);
  err = cudaFree(d_tsc);
  err = cudaFree(d_tss);
  err = cudaFree(d_M);

  free(tcc);
  free(tcs);
  free(tsc);
  free(tss);
  for (int il=0;il<=dim1;il++) {
    free(Tlegq[il]);
    free(Tdleq[il]);
  }
  free(Tlegq);
  free(Tdleq);
  free(h_M);
  // Reset the device and exit
  // cudaDeviceReset causes the driver to clean up all state. While
  // not mandatory in normal operation, it is good practice.  It is also
  // needed to ensure correct operation when the application is being
  // profiled. Calling cudaDeviceReset causes all profile data to be
  // flushed before the application exits
  err = cudaDeviceReset();

  //    return 0;
}




