/*

makematrix.c 

Constructs matrix equation A x = b for local minimization problem 
Exports matrix A and vector b into plain text

Written by Hee Sok Chung at ANL
July 10, 2016

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

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
void cpx(int procid, int **procdec, double **datain, 
         double **cc, double **cs, double **sc, double **ss,
         double ***matrix, 
         double *zeta, double *eta, 
         int **ccidx, int ndatap, int dim1, 
         double zt0, double rr, int *sliceflag);

int main (int argc, char **argv) {
    int nmax, mmax, ndata, dim1, dim2, dim2x, vecdim;
//  Max. number of harmonics in azimuthal direction 
    nmax=200; 
//  Max. number of harmonics in poloidal direction 
    mmax=8; 
//  Matrix and vector dimensions 
    dim1=(nmax+1)*(mmax+1);
    dim2=dim1*dim1;
    dim2x=4*dim1;
    double cterm, rr;
    double **cc, **cs, **sc, **ss;
    double ***matrix, **datain;
    double datatmp;
    double zeta[30], eta[30];
    double zt, et, zt0, bzero;
    double angles[30], radii[30];
    int **ccidx;
    int **mnidx;
    int **procdec;
    int *sinz, *vecidx;
    int i, j, id, idx, md, nd, mdx, ndx, prb, na, ma, qidx;
    int ix, jx, ff, procid, procnum, ndatap, vecc;
    double tp, tx, phi, wgt, iwgt, legq, dleq;
    FILE *input;
    FILE *inputc;
    procnum=11; // number of thread to use 
    printf("Using %d threads\n",procnum);

//  declare vector and matrix elements as pointers 
//  each thread retains its own for maximum speed 
    cc=malloc((procnum+1)*sizeof(double*));
    cs=malloc((procnum+1)*sizeof(double*));
    sc=malloc((procnum+1)*sizeof(double*));
    ss=malloc((procnum+1)*sizeof(double*));
    for (procid=0;procid<=procnum;procid++) {
        cc[procid]=malloc((dim1+2)*sizeof(double));
        cs[procid]=malloc((dim1+2)*sizeof(double));
        sc[procid]=malloc((dim1+2)*sizeof(double));
        ss[procid]=malloc((dim1+2)*sizeof(double));
    }
    matrix=malloc((procnum+1)*sizeof(double**));
    for (procid=0;procid<=procnum;procid++) {
        matrix[procid]=malloc((dim2x+2)*sizeof(double*));
    }
    for (procid=0;procid<=procnum;procid++) {
        for (i=0;i<(dim2x+2);i++) {
            matrix[procid][i]=malloc((dim2x+2)*sizeof(double));
        }
    }
    sinz=malloc((dim2x+2)*sizeof(int)); // for identifying exact zeros 
    vecidx=malloc((dim2x+2)*sizeof(int)); // vector index 
    ccidx=malloc((dim1+2)*sizeof(int*)); // matrix single index -> double index 
    for (i=0;i<(dim1+2);i++) {
        ccidx[i]=malloc(3*sizeof(int));
    }
    mnidx=malloc((dim1+2)*sizeof(int*)); // matrix double index -> single index
    for (i=0;i<(dim1+2);i++) {
        mnidx[i]=malloc((dim1+2)*sizeof(int));
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
    datain=malloc((ndata+1)*sizeof(double*));
    for (id=0;id<=ndata;id++) {
        datain[id]=malloc((25+2)*sizeof(double));
    }
//  distribute azimuthal slices to threads
    ndatap=ndata/procnum+(ndata % procnum != 0);
//    printf("%d\n",procnum);
//    printf("%d\n",ndatap);
    procdec=malloc((procnum+1)*sizeof(int*));
    for (i=0;i<=procnum;i++){ 
        procdec[i]=malloc((ndatap+2)*sizeof(int));
    }
    for (procid=0;procid<=procnum;procid++) {
        for (i=0;i<=ndatap+1;i++) {
            procdec[procid][i]=0;
        }
    }
    procid=1;
    i=1;
    for (id=1;id<=ndata;id++){
        sliceflag[id]=1;
        procdec[procid][i]=id;
        if(procid==procnum) {
            procid=1;
            i++;
        }
        else procid++;
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
    na=0;
    ma=0;
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
            for (procid=0;procid<=procnum;procid++) {
                cc[procid][i]=0.; cs[procid][i]=0.; 
                sc[procid][i]=0.; ss[procid][i]=0.;
            }
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
            for (procid=0;procid<=procnum;procid++) {
                matrix[procid][i][j]=0.;
            }
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

    printf ("\nComputing matrix elements\n");
    omp_set_num_threads(procnum);
#pragma omp parallel
    {
        int idp=omp_get_thread_num();
//  Calculate matrix and vector elements 
    cpx(idp, procdec, datain, cc, cs, sc, ss, matrix, 
         zeta, eta, ccidx, ndatap, dim1, zt0, rr, sliceflag);
    }
    for(i=1;i<=ndata;i++) {
        if (sliceflag[i]!=0) printf("slice %d status = %d\n", i, sliceflag[i]);
    }
//  combine results from individual threads 
    for (i=1;i<=dim2x;i++) {
        for (j=1;j<=dim2x;j++) {
            for (procid=1;procid<=procnum;procid++) {
                matrix[0][i][j]+=matrix[procid][i][j];
            }
        }
    }
    for (i=1;i<=dim1;i++) {
        for (procid=1;procid<=procnum;procid++) {
            cc[0][i]+=cc[procid][i];
            cs[0][i]+=cs[procid][i];
            sc[0][i]+=sc[procid][i];
            ss[0][i]+=ss[procid][i];
        }
    }

    FILE *out0, *out1, *outn, *outd;
    out0 = fopen ("outmatrix.txt", "w"); // matrix elements 
    out1 = fopen ("outvector.txt", "w"); // vector elements 
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
      //fprintf(out1, "%.17g\n", cc[0][i]);
      fwrite(&(cc[0][i]),sizeof(double),1,out1);
      if(sinz[4*(i-1)+2]!=0) {
	vecdim++;
	//fprintf(out1, "%.17g\n", cs[0][i]);
	fwrite(&(cs[0][i]),sizeof(double),1,out1);
      }
      if(sinz[4*(i-1)+3]!=0) {
	vecdim++;
	//fprintf(out1, "%.17g\n", sc[0][i]);
	fwrite(&(sc[0][i]),sizeof(double),1,out1);
      }
      if(sinz[4*(i-1)+4]!=0) {
	vecdim++;
	//fprintf(out1, "%.17g\n", ss[0][i]);
	fwrite(&(ss[0][i]),sizeof(double),1,out1);
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
    //        fprintf(out0, "%.17g\n", matrix[0][vecidx[i]][vecidx[j]]);
	    fwrite(&(matrix[0][vecidx[i]][vecidx[j]]),sizeof(double),1,out0);
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

void cpx(int procid, int **procdec, double **datain, 
         double **cc, double **cs, double **sc, double **ss,
         double ***matrix, 
         double *zeta, double *eta, 
         int **ccidx, int ndatap, int dim1, 
         double zt0, double rr, int *sliceflag) {
    double datatmp, zt, et, wgt, phi;
    double *tcc, *tcs, *tsc, *tss;
    tcc=(double *)malloc(sizeof*tcc*(dim1+2));
    tcs=(double *)malloc(sizeof*tcs*(dim1+2));
    tsc=(double *)malloc(sizeof*tsc*(dim1+2));
    tss=(double *)malloc(sizeof*tss*(dim1+2));
    int md, nd, i, id, prb, j, ix, jx, ixx;
    double legq, dleq;
    printf("Starting thread %d\n",procid);

    //Calculate legq and dleq table
    double **Tlegq;
    double **Tdleq;
    Tlegq=malloc((dim1+1)*sizeof(double*));
    Tdleq=malloc((dim1+1)*sizeof(double*));
    for (int il=0;il<=dim1;il++) {
      Tlegq[il]=malloc(26*sizeof(double));
      Tdleq[il]=malloc(26*sizeof(double));
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

    for (ixx=1;ixx<=ndatap;ixx++) { // loop over azimuthal slices 
        id=procdec[procid+1][ixx];
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
                tcc[i]=-Sinh_zt*Sin_et/rr*cos(nd*phi)*cos(md*et)*  \
                    (Sinh_zt/2./wgt*legq+wgt*dleq) \
                    -(1.-Cosh_zt*Cos_et)/rr* \
                    (-md*cos(nd*phi)*sin(md*et)*wgt*legq \
                     +cos(nd*phi)*cos(md*et)*Sin_et/2./wgt*legq);
// CS(m,n) coefficient 
                tcs[i]=-Sinh_zt*Sin_et/rr*sin(nd*phi)*cos(md*et)* \
                    (Sinh_zt/2./wgt*legq+wgt*dleq) \
                    -(1.-Cosh_zt*Cos_et)/rr* \
                    (-md*sin(nd*phi)*sin(md*et)*wgt*legq \
                     +sin(nd*phi)*cos(md*et)*Sin_et/2./wgt*legq);
// SC(m,n) coefficient 
                tsc[i]=-Sinh_zt*Sin_et/rr*cos(nd*phi)*sin(md*et)* \
                    (Sinh_zt/2./wgt*legq+wgt*dleq) \
                    -(1.-Cosh_zt*Cos_et)/rr* \
                    (md*cos(nd*phi)*cos(md*et)*wgt*legq \
                     +cos(nd*phi)*sin(md*et)*Sin_et/2./wgt*legq);
// SS(m,n) coefficient 
                tss[i]=-Sinh_zt*Sin_et/rr*sin(nd*phi)*sin(md*et)* \
                    (Sinh_zt/2./wgt*legq+wgt*dleq) \
                    -(1.-Cosh_zt*Cos_et)/rr* \
                    (md*sin(nd*phi)*cos(md*et)*wgt*legq \
                     +sin(nd*phi)*sin(md*et)*Sin_et/2./wgt*legq);
            } // define tcc loop end
            for (i=1;i<=dim1;i++){
                for (j=1;j<=dim1;j++) {
// Double indices 
                    ix=i-1; jx=j-1;
                    matrix[procid][4*ix+1][4*jx+1]+=tcc[i]*tcc[j]; 
                    matrix[procid][4*ix+2][4*jx+1]+=tcs[i]*tcc[j]; 
                    matrix[procid][4*ix+3][4*jx+1]+=tsc[i]*tcc[j]; 
                    matrix[procid][4*ix+4][4*jx+1]+=tss[i]*tcc[j]; 

                    matrix[procid][4*ix+1][4*jx+2]+=tcc[i]*tcs[j]; 
                    matrix[procid][4*ix+2][4*jx+2]+=tcs[i]*tcs[j]; 
                    matrix[procid][4*ix+3][4*jx+2]+=tsc[i]*tcs[j]; 
                    matrix[procid][4*ix+4][4*jx+2]+=tss[i]*tcs[j]; 

                    matrix[procid][4*ix+1][4*jx+3]+=tcc[i]*tsc[j]; 
                    matrix[procid][4*ix+2][4*jx+3]+=tcs[i]*tsc[j]; 
                    matrix[procid][4*ix+3][4*jx+3]+=tsc[i]*tsc[j]; 
                    matrix[procid][4*ix+4][4*jx+3]+=tss[i]*tsc[j]; 

                    matrix[procid][4*ix+1][4*jx+4]+=tcc[i]*tss[j]; 
                    matrix[procid][4*ix+2][4*jx+4]+=tcs[i]*tss[j]; 
                    matrix[procid][4*ix+3][4*jx+4]+=tsc[i]*tss[j]; 
                    matrix[procid][4*ix+4][4*jx+4]+=tss[i]*tss[j]; 
// Quadratic terms 
                }
            } // define matrix loop end
// Linear terms 
            for (i=1;i<=dim1;i++){
                cc[procid][i]+=tcc[i]*datatmp;
                cs[procid][i]+=tcs[i]*datatmp;
                sc[procid][i]+=tsc[i]*datatmp;
                ss[procid][i]+=tss[i]*datatmp;
            } // define vector loop end 
// for each data end
        } // probe loop end 
//    printf ("\nSlice %d (phi=%lf) done\n", id, phi);
    printf ("\nThread %d, localid %d of %d done\n", procid, ixx, ndatap);
    }} // id end 
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

//    return 0;
}




