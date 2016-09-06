/*

bfield.c 

Computes B fields in toroidal coordinates with given coefficients 
and calculates rms deviation from data

Written by Hee Sok Chung at ANL
July 10, 2016

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <TTree.h>
//#include <TFile.h>
#include <iostream>
#include <cuda_runtime.h>

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

typedef struct B_Struct{
  double Prb[25];
}B_Struct;

using namespace std;

__global__ void thread_call(const double *cc,const double *cs,const double *sc,const double *ss,const double *zeta,const double *eta,const double *Tlegq,const double *Tdleq,const int *midx,const int *nidx,const int*ccidx,const double *dataread,double *dataoutZ,double *dataoutR,double *dataoutPhi,const double rr,const double bzero,const int dim1,int nphi);

int main () {
  double *cc, *cs, *sc, *ss, *lsigma;
  int *ccidx;
  int *midx;
  int *nidx;
  double *Tlegq;
  double *Tdleq;
  double *zeta, *eta;
  double * dataread;
  double * dataoutZ;
  double * dataoutR;
  double * dataoutPhi;
  cudaError_t err = cudaSuccess;

  zeta = (double*)malloc(30*sizeof(double));
  eta = (double*)malloc(30*sizeof(double));

    int nmax, mmax, nphi, dim1, dim2;
    FILE *input, *iname, *datafile, *datac , *output;
    int Nthread = 8;
//    input = fopen ("fitcoeffN.txt", "r");
    input = fopen ("fitc.txt", "r");
    iname = fopen ("outnames.txt", "r");
    output = fopen ("output52full.txt","w");
    fscanf (iname, "%d", &mmax); // read m Maximum
    fscanf (iname, "%d", &nmax); // read n Maximum
//    fscanf (input, "%d", &mmax); // read m Maximum
//    fscanf (input, "%d", &nmax); // read n Maximum
    printf("MMax %d, NMax %d \n",mmax, nmax);
    dim1=(nmax+1)*(mmax+1); // coefficients dimensions 
    double cterm, rr;
//    int ccidx[mmax+2][nmax+2];
    double zt, et, zt0, bzero, br, bz, bphi, bsize;
    double difrms, dzrms, ffrms, rmsd;
    double angles[30], radii[30], bfield[30][5];
    int i, j, id, idx, md, nd, mdx, ndx, prb, na, ma, qidx, vid;
    double datatmp, tp, tx, phi, wgt, iwgt, legq, dleq;
    double phidz, phidp, phide, phide1, phide2;
    B_Struct B_Measure;
    B_Struct B_Fit;
    B_Struct B_ZFit;
    B_Struct B_RFit;
    B_Struct B_PhiFit;
    cc=(double *)malloc(sizeof*cc*(dim1+2));
    cs=(double *)malloc(sizeof*cs*(dim1+2));
    sc=(double *)malloc(sizeof*sc*(dim1+2));
    ss=(double *)malloc(sizeof*ss*(dim1+2));
    lsigma=(double *)malloc(sizeof*lsigma*(nmax+2));
    int ccidx_dim = dim1+2;
    ccidx=(int *)malloc((dim1+2)*(dim1+2)*sizeof(int));
    if(ccidx==NULL){
        printf("out of memory\n");
    }
    midx=(int *)malloc((dim1+2)*sizeof(int*));
    nidx=(int *)malloc((dim1+2)*sizeof(int*));


    bzero=61.789; // average B-field 
    zt0=6.5;          // eigenfunction normalization point
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
    for (i=1;i<=25;i++) {
        printf("probe %d location r=%.17g, theta=%.17g \n",i,radii[i],angles[i]);
    }

//  calculate toroidal coordinates 
    rr=7111.5; // toroid center 
    printf ("\nComputing coordinates\n");
    for (i=1;i<=25;i++) {
        zeta[i]=zetaf(radii[i]*sin(angles[i])+7112., \
                      -radii[i]*cos(angles[i]),rr);
        eta[i]=etaf(radii[i]*sin(angles[i])+7112., \
                    -radii[i]*cos(angles[i]),rr);
    }
    /* initialize arrays */
    na=0;
    ma=0;
    i=0;
    nd=0;
    md=0;
    lsigma[0]=1.;
    lsigma[nmax]=0.;
    for (nd=0;nd<=nmax;nd++) {
        for (md=0;md<=mmax;md++) {
            i+=1;
            ccidx[md*ccidx_dim+nd]=i;
            midx[i]=md;
            nidx[i]=nd;
            cc[i]=0.; cs[i]=0.; sc[i]=0.; ss[i]=0.;
        }
        if(nd>0&&nd<nmax) {
            lsigma[nd]=sin(M_PI*(nd*1.)/(nmax*1.));
            lsigma[nd]=lsigma[nd]/(M_PI*(nd*1.)/(nmax*1.));
        }

    }
    for (idx=1;idx<=dim1*4;idx++) { // read coefficients
        fscanf (iname, "%d \t %d \t %d", &id, &ma, &na); 
//        fscanf (input, "%d \t %d \t %d \t %lg", &id, &ma, &na, &tp); 
        fscanf (input, "%lg", &tp); 
//        printf("%d, %d, %lg \n", id, ccidx[ma*ccidx_dim+na], tp);
//        if (tp!=tp) printf("%d, %d, %lg \n", id, ccidx[ma][na], tp);
        if (id==1) {
//            printf("cc \n");
            cc[ccidx[ma*ccidx_dim+na]]=tp;
        }
        else if (id==2) {
//            printf("cs \n");
            cs[ccidx[ma*ccidx_dim+na]]=tp;
        }
        else if (id==3) {
//            printf("sc \n");
            sc[ccidx[ma*ccidx_dim+na]]=tp;
        }
        else if (id==4) {
//            printf("ss \n");
            ss[ccidx[ma*ccidx_dim+na]]=tp;
        }
        else {
            printf("error \n");
        }
//        printf("%.20lg \n",tp);
    }
    fclose(input); // close data 
    fclose(iname); // close data 
    printf ("\nComputing b-fields\n");
    nphi=0;
    datac=fopen("data52.txt", "r"); // open data for counting
    while(fscanf(datac, "%lg", &datatmp)>0) {
    nphi++;
    }
    nphi=nphi/26; // azimuthal angle + 25 probes 
    fclose(datac); // close data file used for counting 

    //Loading data first
    dataread = (double *)malloc(sizeof(double)*26*nphi);
    dataoutZ = (double *)malloc(sizeof(double)*26*nphi);
    dataoutR = (double *)malloc(sizeof(double)*26*nphi);
    dataoutPhi = (double *)malloc(sizeof(double)*26*nphi);
    double tmp;
    datafile = fopen ("data52.txt", "r");
    for (id=0;id<nphi;id++) { // loop over azimuthal slices 
      fscanf (datafile, "%lg", &tmp); // read azimuthal angle in degrees 
      dataread[id*26]=tmp*M_PI/180.;
      for (prb=1;prb<=25;prb++){ // loop over 25 probes 
	fscanf (datafile, "%lg", &tmp); 
	dataread[id*26+prb]=tmp*0.001+61.7400000;
      }
    }
    fclose(datafile);
    //Zero output
    for (id=0;id<nphi;id++) { // loop over azimuthal slices 
      for (prb=1;prb<=25;prb++){ // loop over 25 probes 
	dataoutZ[id*26+prb]=0.0;
	dataoutR[id*26+prb]=0.0;
	dataoutPhi[id*26+prb]=0.0;
      }
    }
    //Calculate legq and dleq table
    Tlegq=(double *)malloc((dim1+1)*26*sizeof(double));
    Tdleq=(double *)malloc((dim1+1)*26*sizeof(double));

    for (prb=0;prb<=25;prb++){ // loop over 25 probes 
      for (i=0;i<=dim1;i++){
	Tlegq[i*26+prb]=0.0;
	Tdleq[i*26+prb]=0.0;
      }
    }
    for (prb=1;prb<=25;prb++){ // loop over 25 probes 
      for (i=1;i<=dim1;i++){
	int md=midx[i]; // m index 
	int nd=nidx[i]; // n index 
	/* LegendreQ at probe */
	Tlegq[i*26+prb]=LegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
	/* Derivative of LegendreQ at probe */
	Tdleq[i*26+prb]=DLegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
      }
    }

    //
    difrms=0.;
    dzrms=0.;
    ffrms=0.;

    //Allocate memory in device
    //Function Talbe
    double *d_legq = NULL;
    double *d_dleq = NULL;
    int f_size = sizeof(double)*26*(dim1+1);
    err = cudaMalloc((void **)&d_legq, f_size);
    err = cudaMalloc((void **)&d_dleq, f_size);

    //Cooridinates
    double *d_zeta = NULL;
    double *d_eta = NULL;
    int c_size = sizeof(double)*30;
    err = cudaMalloc((void **)&d_zeta, c_size);
    err = cudaMalloc((void **)&d_eta, c_size);

    //data storage;
    double *d_Data = NULL;
    double *d_DataoutZ = NULL;
    double *d_DataoutR = NULL;
    double *d_DataoutPhi = NULL;
    size_t sizeData = sizeof(double)*26*nphi;
    err = cudaMalloc((void **)&d_Data, sizeData);
    err = cudaMalloc((void **)&d_DataoutZ, sizeData);
    err = cudaMalloc((void **)&d_DataoutR, sizeData);
    err = cudaMalloc((void **)&d_DataoutPhi, sizeData);

    //ccidx
    int *d_ccidx = NULL;
    int size_ccidx = (dim1+2)*(dim1+2)*sizeof(int);
    err = cudaMalloc((void **)&d_ccidx, size_ccidx);

    int *d_midx = NULL;
    int size_midx = (dim1+2)*sizeof(int);
    err = cudaMalloc((void **)&d_midx, size_midx);

    int *d_nidx = NULL;
    int size_nidx = (dim1+2)*sizeof(int);
    err = cudaMalloc((void **)&d_nidx, size_nidx);

    //Vectors
    double *d_cc = NULL;
    double *d_cs = NULL;
    double *d_sc = NULL;
    double *d_ss = NULL;
    int d_size = sizeof(double)*(dim1+2);
    err = cudaMalloc((void **)&d_cc, d_size);
    err = cudaMalloc((void **)&d_cs, d_size);
    err = cudaMalloc((void **)&d_sc, d_size);
    err = cudaMalloc((void **)&d_ss, d_size);

    if (err != cudaSuccess)
    {
      fprintf(stderr, "Failed to allocate device memory (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }


    //copy matrix memory to device
    err = cudaMemcpy(d_legq, Tlegq, f_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_dleq, Tdleq, f_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_zeta, zeta, c_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_eta, eta, c_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_Data, dataread, sizeData, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_DataoutZ, dataoutZ, sizeData, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_DataoutR, dataoutR, sizeData, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_DataoutPhi, dataoutPhi, sizeData, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_ccidx, ccidx, size_ccidx, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_midx, midx, size_midx, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_nidx, nidx, size_nidx, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_cc, cc, d_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_cs, cs, d_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_sc, sc, d_size, cudaMemcpyHostToDevice);
    err = cudaMemcpy(d_ss, ss, d_size, cudaMemcpyHostToDevice);



    //Start parallel computing
    dim3 DimBlock (16,16);
    dim3 DimGrid (nphi/16+1, 26/16+1);
    printf("CUDA kernel launch with %d blocks of %d threads\n", DimGrid.x * DimGrid.y, 256);

    thread_call<<<DimGrid, DimBlock>>>(d_cc,d_cs,d_sc,d_ss,d_zeta,d_eta,d_legq,d_dleq,d_midx,d_nidx,d_ccidx,d_Data,d_DataoutZ,d_DataoutR,d_DataoutPhi,rr,bzero,dim1,nphi);
    err = cudaGetLastError();

    if (err != cudaSuccess)
    {
      fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
      exit(EXIT_FAILURE);
    }
    printf("GPU jobs done.\n");


    //End parallel computing
    err = cudaMemcpy(dataoutZ, d_DataoutZ, sizeData, cudaMemcpyDeviceToHost);
    err = cudaMemcpy(dataoutR, d_DataoutR, sizeData, cudaMemcpyDeviceToHost);
    err = cudaMemcpy(dataoutPhi, d_DataoutPhi, sizeData, cudaMemcpyDeviceToHost);
    
    //open Root tree
/*    TFile *outfile = new TFile("RootOut52.root","recreate");
    TTree *Tree_Measured = new TTree ("Tree_Measured", "Measured field");
    Tree_Measured->Branch("Phi",&phi,"Phi/D");
    Tree_Measured->Branch("BField",&B_Measure,"Prb1/D:Prb2:Prb3:Prb4:Prb5:Prb6:Prb7:Prb8:Prb9:Prb10:Prb11:Prb12:Prb13:Prb14:Prb15:Prb16:Prb17:Prb18:Prb19:Prb20:Prb21:Prb22:Prb23:Prb24:Prb25");
    TTree *Tree_Fit = new TTree ("Tree_Fit", "Fitted field");
    Tree_Fit->Branch("Phi",&phi,"Phi/D");
    Tree_Fit->Branch("BField",&B_Fit,"Prb1/D:Prb2:Prb3:Prb4:Prb5:Prb6:Prb7:Prb8:Prb9:Prb10:Prb11:Prb12:Prb13:Prb14:Prb15:Prb16:Prb17:Prb18:Prb19:Prb20:Prb21:Prb22:Prb23:Prb24:Prb25");
    Tree_Fit->Branch("BFieldZ",&B_Fit,"Prb1/D:Prb2:Prb3:Prb4:Prb5:Prb6:Prb7:Prb8:Prb9:Prb10:Prb11:Prb12:Prb13:Prb14:Prb15:Prb16:Prb17:Prb18:Prb19:Prb20:Prb21:Prb22:Prb23:Prb24:Prb25");
    Tree_Fit->Branch("BFieldR",&B_Fit,"Prb1/D:Prb2:Prb3:Prb4:Prb5:Prb6:Prb7:Prb8:Prb9:Prb10:Prb11:Prb12:Prb13:Prb14:Prb15:Prb16:Prb17:Prb18:Prb19:Prb20:Prb21:Prb22:Prb23:Prb24:Prb25");
    Tree_Fit->Branch("BFieldPhi",&B_Fit,"Prb1/D:Prb2:Prb3:Prb4:Prb5:Prb6:Prb7:Prb8:Prb9:Prb10:Prb11:Prb12:Prb13:Prb14:Prb15:Prb16:Prb17:Prb18:Prb19:Prb20:Prb21:Prb22:Prb23:Prb24:Prb25");
    */
    for (id=0;id<nphi;id++) { // loop over azimuthal slices 
      phi=dataread[id*26];
      for (prb=1;prb<=25;prb++){ // loop over 25 probes 
	bz=dataoutZ[id*26+prb];
	br=dataoutR[id*26+prb];
	bphi=dataoutPhi[id*26+prb];
	bsize=sqrt(br*br+bz*bz+bphi*bphi);
	//            printf("%Lf \n",bsize);
	bfield[prb][4]=bsize;
	bfield[prb][1]=bz;
	bfield[prb][2]=br;
	bfield[prb][3]=bphi;
	B_Fit.Prb[prb-1]=bsize;
	B_ZFit.Prb[prb-1]=bz;
	B_RFit.Prb[prb-1]=br;
	B_PhiFit.Prb[prb-1]=bphi;
	B_Measure.Prb[prb-1]=dataread[id*26+prb];
	//Output
	fprintf(output,"%.17g %d %.17g %.17g %.17g %.17g %.17g\n",phi,prb,dataread[id*26+prb],bfield[prb][1],bfield[prb][2],bfield[prb][3],bfield[prb][4]);
	rmsd=(dataread[id*26+prb]-bfield[prb][4])*(dataread[id*26+prb]-bfield[prb][4]);
	difrms+=rmsd;
	rmsd=rmsd/bzero/bzero;
	dzrms+=(dataread[id*26+prb]-bz)*(dataread[id*26+prb]-bz);
	ffrms+=(dataread[id*26+prb]-bzero)*(dataread[id*26+prb]-bzero);
      }
     // Tree_Measured->Fill();
     // Tree_Fit->Fill();
    }
    //Tree_Measured->Write();
    //Tree_Fit->Write();
    //outfile->Close();

    difrms=difrms/bzero/bzero/nphi/25;
    dzrms=dzrms/bzero/bzero/nphi/25;
    ffrms=ffrms/bzero/bzero/nphi/25;
    difrms=sqrt(difrms)*1000000.;
    dzrms=sqrt(dzrms)*1000000.;
    ffrms=sqrt(ffrms)*1000000.;
    printf("RMS fluctuation = %.17g ppm \n", ffrms);
    printf("RMS difference (lin. approx.) = %.17g ppm \n", dzrms);
    printf("RMS difference (real) = %.17g ppm \n", difrms);

    free(cc);
    free(cs);
    free(sc);
    free(ss);
    free(ccidx);
    free(midx);
    free(nidx);
    free(Tlegq);
    free(Tdleq);
    free(dataread);
    free(dataoutZ);
    free(dataoutR);
    free(dataoutPhi);

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
		f1f2=1./err; // to return inf or nan
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
//    lq=pow(M_PI,.5)*exp(lgammal(m+n+.5))/(pow(2.,m+.5)*exp(lgammal(m+1.)));
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
    zetax=atanhl(2.*rho*r0/(rho*rho+r0*r0+z*z));
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


//Function for each thread call
__global__ void thread_call(const double *cc,const double *cs,const double *sc,const double *ss,const double *zeta,const double *eta,const double *Tlegq,const double *Tdleq,const int *midx,const int *nidx,const int*ccidx,const double *dataread,double *dataoutZ,double *dataoutR,double *dataoutPhi,const double rr,const double bzero,const int dim1,int nphi){
  //Get index
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  int prb = blockDim.y * blockIdx.y + threadIdx.y;
  int ccidx_dim = dim1+2;
  if(id<nphi && prb<26 && prb>0){
    double phi=dataread[id*26];
    double zt=zeta[prb];  // zeta coordinate at probe 
    double et=eta[prb];  // eta coordinate at probe 
    //            printf("%Lf, %Lf \n", zt, et);
    double wgt=sqrt(cosh(zt)-cos(et)); // weight func in toroidal coordinates 
    double br=0.; 
    double bz=bzero; 
    double bphi=0.;
    for (int i=1;i<=dim1;i++){
      int md=midx[i]; // m index 
      int nd=nidx[i]; // n index 
      //                printf("%d, %d \n", md, nd);
      /* LegendreQ at probe */
      double legq=Tlegq[i*26+prb];
      /* Derivative of LegendreQ at probe */
      double dleq=Tdleq[i*26+prb];
      //                printf("(%d, %d, %Lf) : %Lf, %Lf \n", md,nd,zt,legq, dleq);
      double phidz=cc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*cos(md*et)+  \
		   sc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*sin(md*et)+  \
		   cs[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*cos(md*et)+  \
		   ss[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*sin(md*et);
      phidz=phidz*(sinh(zt)/2./wgt*legq + wgt*dleq);
      double phidp=-cc[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*cos(md*et)-  \
		   sc[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*sin(md*et)+  \
		   cs[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*cos(md*et)+  \
		   ss[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*sin(md*et);
      phidp=phidp*nd*wgt*legq;
      double phide1=cc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*cos(md*et)+  \
		    sc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*sin(md*et)+  \
		    cs[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*cos(md*et)+  \
		    ss[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*sin(md*et);
      phide1=phide1*sin(et)/2./wgt*legq;
      double phide2=-cc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*sin(md*et)+  \
		    sc[ccidx[md*ccidx_dim+nd]]*cos(nd*phi)*cos(md*et)-  \
		    cs[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*sin(md*et)+  \
		    ss[ccidx[md*ccidx_dim+nd]]*sin(nd*phi)*cos(md*et);
      phide2=phide2*md*wgt*legq;
      double phide=phide1+phide2;
      br=br+sinh(zt)/rr*((1.-cosh(zt)*cos(et))/sinh(zt)*phidz-sin(et)*phide);
      bz=bz+sinh(zt)/rr*(-sin(et)*phidz-(1.-cosh(zt)*cos(et))/sinh(zt)*phide);
      bphi=bphi+(cosh(zt)-cos(et))/(rr*sinh(zt))*phidp;
    }
    dataoutZ[id*26+prb]=bz;
    dataoutR[id*26+prb]=br;
    dataoutPhi[id*26+prb]=bphi;
  }
}







