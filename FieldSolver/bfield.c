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

/* hypergeometric function 1F2 */
long double h1f2(long double a, long double b, long double c, long double z);

/* LegendreQ(m-1/2,n,cosh(x)) */
long double LegendreQ(long double m, long double n, long double z);

/* derivative of LegendreQ(m-1/2,n,cosh(x)) */
long double DLegendreQ(int m, int n, long double z);

/* coordinate transformation r, z -> zeta */
long double zetaf(long double rho, long double z, long double r0);

/* coordinate transformation r, z -> eta */
long double etaf(long double rho, long double z, long double r0);

int main () {
    int nmax, mmax, nphi, dim1, dim2;
    FILE *input, *iname, *datafile, *datac;
//    input = fopen ("fitcoeffN.txt", "r");
    input = fopen ("fitc.txt", "r");
    iname = fopen ("outnames.txt", "r");
//    datafile = fopen ("data40.txt", "r");
    datafile = fopen ("data52.txt", "r");
    fscanf (iname, "%d", &mmax); // read m Maximum
    fscanf (iname, "%d", &nmax); // read n Maximum
//    fscanf (input, "%d", &mmax); // read m Maximum
//    fscanf (input, "%d", &nmax); // read n Maximum
    printf("MMax %d, NMax %d \n",mmax, nmax);
    dim1=(nmax+1)*(mmax+1); // coefficients dimensions 
    long double cterm, rr;
    long double *cc, *cs, *sc, *ss, *lsigma;
    int **ccidx;
    int *midx;
    int *nidx;
//    int ccidx[mmax+2][nmax+2];
    long double zeta[30], eta[30];
    long double zt, et, zt0, bzero, br, bz, bphi, bsize;
    long double dataread[30], difrms, dzrms, ffrms, rmsd;
    long double angles[30], radii[30], bfield[30][5];
    int i, j, id, idx, md, nd, mdx, ndx, prb, na, ma, qidx, vid;
    long double datatmp, tp, tx, phi, wgt, iwgt, legq, dleq;
    long double phidz, phidp, phide, phide1, phide2;
    cc=(long double *)malloc(sizeof*cc*(dim1+2));
    cs=(long double *)malloc(sizeof*cs*(dim1+2));
    sc=(long double *)malloc(sizeof*sc*(dim1+2));
    ss=(long double *)malloc(sizeof*ss*(dim1+2));
    lsigma=(long double *)malloc(sizeof*lsigma*(nmax+2));
    ccidx=malloc((dim1+2)*sizeof(int*));
    if(ccidx==NULL){
        printf("out of memory\n");
    }
    for (i=0;i<(dim1+2);i++) {
        ccidx[i]=malloc((dim1+2)*sizeof(int));
        if (ccidx[i]==NULL) {
            printf("out of memory\n");
        }
    }
    midx=malloc((dim1+2)*sizeof(int*));
    nidx=malloc((dim1+2)*sizeof(int*));


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
        printf("probe %d location r=%Lf, theta=%Lf \n",i,radii[i],angles[i]);
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
            ccidx[md][nd]=i;
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
//        fscanf (input, "%d \t %d \t %d \t %Lf", &id, &ma, &na, &tp); 
        fscanf (input, "%Lf", &tp); 
        printf("%d, %d, %Lf \n", id, ccidx[ma][na], tp);
//        if (tp!=tp) printf("%d, %d, %Lf \n", id, ccidx[ma][na], tp);
        if (id==1) {
//            printf("cc \n");
            cc[ccidx[ma][na]]=tp;
        }
        else if (id==2) {
//            printf("cs \n");
            cs[ccidx[ma][na]]=tp;
        }
        else if (id==3) {
//            printf("sc \n");
            sc[ccidx[ma][na]]=tp;
        }
        else if (id==4) {
//            printf("ss \n");
            ss[ccidx[ma][na]]=tp;
        }
        else {
            printf("error \n");
        }
//        printf("%.20Lf \n",tp);
    }
    fclose(input); // close data 
    fclose(iname); // close data 
    printf ("\nComputing b-fields\n");
    nphi=0;
//    datac=fopen("data40.txt", "r"); // open data for counting
    datac=fopen("data52.txt", "r"); // open data for counting
    while(fscanf(datac, "%Lf", &datatmp)>0) {
    nphi++;
    }
    nphi=nphi/26; // azimuthal angle + 25 probes 
    fclose(datac); // close data file used for counting 
    difrms=0.;
    dzrms=0.;
    ffrms=0.;
    for (id=1;id<=nphi;id++) { // loop over azimuthal slices 
        fscanf (datafile, "%Lf", &dataread[0]); // read azimuthal angle in degrees 
        phi=dataread[0]*M_PI/180.;
        for (prb=1;prb<=25;prb++){ // loop over 25 probes 
            fscanf (datafile, "%Lf", &dataread[prb]); 
            dataread[prb]=dataread[prb]*0.001+61.7400000;
            zt=zeta[prb];  // zeta coordinate at probe 
            et=eta[prb];  // eta coordinate at probe 
//            printf("%Lf, %Lf \n", zt, et);
            wgt=sqrt(cosh(zt)-cos(et)); // weight func in toroidal coordinates 
            br=0.; bz=bzero; bphi=0.;
            for (i=1;i<=dim1;i++){
                md=midx[i]; // m index 
                nd=nidx[i]; // n index 
//                printf("%d, %d \n", md, nd);
// LegendreQ at probe
                legq=LegendreQ(md,nd,zt)/LegendreQ(md,nd,zt0);
// Derivative of LegendreQ at probe
                dleq=DLegendreQ(md,nd,zt)/LegendreQ(md,nd,zt0);
//                printf("(%d, %d, %Lf) : %Lf, %Lf \n", md,nd,zt,legq, dleq);
                phidz=cc[ccidx[md][nd]]*cos(nd*phi)*cos(md*et)+  \
                      sc[ccidx[md][nd]]*cos(nd*phi)*sin(md*et)+  \
                      cs[ccidx[md][nd]]*sin(nd*phi)*cos(md*et)+  \
                      ss[ccidx[md][nd]]*sin(nd*phi)*sin(md*et);
                phidz=phidz*(sinh(zt)/2./wgt*legq + wgt*dleq);
                phidp=-cc[ccidx[md][nd]]*sin(nd*phi)*cos(md*et)-  \
                      sc[ccidx[md][nd]]*sin(nd*phi)*sin(md*et)+  \
                      cs[ccidx[md][nd]]*cos(nd*phi)*cos(md*et)+  \
                      ss[ccidx[md][nd]]*cos(nd*phi)*sin(md*et);
                phidp=phidp*nd*wgt*legq;
                phide1=cc[ccidx[md][nd]]*cos(nd*phi)*cos(md*et)+  \
                       sc[ccidx[md][nd]]*cos(nd*phi)*sin(md*et)+  \
                       cs[ccidx[md][nd]]*sin(nd*phi)*cos(md*et)+  \
                       ss[ccidx[md][nd]]*sin(nd*phi)*sin(md*et);
                phide1=phide1*sin(et)/2./wgt*legq;
                phide2=-cc[ccidx[md][nd]]*cos(nd*phi)*sin(md*et)+  \
                       sc[ccidx[md][nd]]*cos(nd*phi)*cos(md*et)-  \
                       cs[ccidx[md][nd]]*sin(nd*phi)*sin(md*et)+  \
                       ss[ccidx[md][nd]]*sin(nd*phi)*cos(md*et);
                phide2=phide2*md*wgt*legq;
                phide=phide1+phide2;
                br=br+sinh(zt)/rr*((1.-cosh(zt)*cos(et))/sinh(zt)*phidz-sin(et)*phide);
                bz=bz+sinh(zt)/rr*(-sin(et)*phidz-(1.-cosh(zt)*cos(et))/sinh(zt)*phide);
                bphi=bphi+(cosh(zt)-cos(et))/(rr*sinh(zt))*phidp;
            }
            bsize=sqrt(br*br+bz*bz+bphi*bphi);
//            printf("%Lf \n",bsize);
            bfield[prb][4]=bsize;
            bfield[prb][1]=bz;
            bfield[prb][2]=br;
            bfield[prb][3]=bphi;
            rmsd=(dataread[prb]-bfield[prb][4])*(dataread[prb]-bfield[prb][4]);
            difrms+=rmsd;
            rmsd=rmsd/bzero/bzero;
            dzrms+=(dataread[prb]-bz)*(dataread[prb]-bz);
            ffrms+=(dataread[prb]-bzero)*(dataread[prb]-bzero);
        }
        printf("Slice %d (phi=%Lf) of %d done (%Lf, %Lf, %Lf) \n",
                id, phi, nphi, difrms, dzrms, ffrms);
    }
    difrms=difrms/bzero/bzero/nphi/25;
    dzrms=dzrms/bzero/bzero/nphi/25;
    ffrms=ffrms/bzero/bzero/nphi/25;
    difrms=sqrt(difrms)*1000000.;
    dzrms=sqrt(dzrms)*1000000.;
    ffrms=sqrt(ffrms)*1000000.;
    printf("RMS fluctuation = %Lf ppm \n", ffrms);
    printf("RMS difference (lin. approx.) = %Lf ppm \n", dzrms);
    printf("RMS difference (real) = %Lf ppm \n", difrms);
    /* calculate coeffs */
/*
    cterm=0.;
    printf ("\nComputing b-fields\n");
    for (id=1;id<=ndata;id++) { // loop over azimuthal slices 
        fscanf (input, "%Lf", &datatmp); // read azimuthal angle in degrees 
        phi=datatmp*M_PI/180.;
        for (prb=1;prb<=25;prb++){ // loop over 25 probes 
            zt=zeta[prb];  // zeta coordinate at probe 
            et=eta[prb];  // eta coordinate at probe 
            wgt=sqrt(cosh(zt)-cos(et)); // weight func in toroidal coordinates 
            br=0.; bz=0.; bphi=0.;
            for (i=1;i<=dim1;i++){
                md=ccidx[i][1]; // m index 
                nd=ccidx[i][2]; // n index 
// LegendreQ at probe
                legq=LegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
// Derivative of LegendreQ at probe
                dleq=DLegendreQ(md,nd,zeta[prb])/LegendreQ(md,nd,zt0);
            }
            idx=0;
// for each data end
        }
    printf ("\nSlice %d (phi=%Lf) done\n", id, phi);
    }
*/
/*

    FILE *out0;

    out0 = fopen ("chisqc.m", "w");

    fprintf(out0, "NMax=%d;\n", nmax);
    fprintf(out0, "MMax=%d;\n", mmax);
    fprintf(out0, "csqc0=%.20Lf;\n", cterm);
    fprintf(out0, "csqc=(%.20Lf)+\n", cterm);

    for (i=1;i<=dim1;i++) {
        fprintf(out0, "(%.20Lf)*CC[%d,%d]+ \n", \
                cc[i],ccidx[i][1],ccidx[i][2]);
        if (ccidx[i][2]!=0) fprintf(out0, "(%.20Lf)*CS[%d,%d]+ \n", \
                cs[i],ccidx[i][1],ccidx[i][2]);
        if (ccidx[i][1]!=0) fprintf(out0, "(%.20Lf)*SC[%d,%d]+ \n", \
                sc[i],ccidx[i][1],ccidx[i][2]);
        if (ccidx[i][1]!=0&&ccidx[i][2]!=0) fprintf(out0, \
                "(%.20Lf)*SS[%d,%d]+ \n", \
                ss[i],ccidx[i][1],ccidx[i][2]);
    }
    for (i=1;i<=dim2;i++) {
        fprintf(out0, "(%.20Lf)*CC[%d,%d]*CC[%d,%d]+ \n",\
        cccc[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*CC[%d,%d]*CS[%d,%d]+ \n",\
        cccs[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][3]!=0) fprintf(out0, \
                    "(%.20Lf)*CC[%d,%d]*SC[%d,%d]+ \n",\
        ccsc[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][3]!=0&&ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*CC[%d,%d]*SS[%d,%d]+ \n",\
        ccss[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][2]!=0&&ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*CS[%d,%d]*CS[%d,%d]+ \n",\
        cscs[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][2]!=0&&ccccidx[i][3]!=0) fprintf(out0, \
                    "(%.20Lf)*CS[%d,%d]*SC[%d,%d]+ \n",\
        cssc[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][2]!=0&&ccccidx[i][3]!=0&&ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*CS[%d,%d]*SS[%d,%d]+ \n",\
        csss[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][1]!=0&&ccccidx[i][3]!=0) fprintf(out0, \
                    "(%.20Lf)*SC[%d,%d]*SC[%d,%d]+ \n",\
        scsc[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][1]!=0&&ccccidx[i][3]!=0&&ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*SC[%d,%d]*SS[%d,%d]+ \n",\
        scss[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][1]!=0&&ccccidx[i][2]&& \
            ccccidx[i][3]!=0&&ccccidx[i][4]!=0) fprintf(out0, \
                    "(%.20Lf)*SS[%d,%d]*SS[%d,%d]",\
        ssss[i], ccccidx[i][1],ccccidx[i][2],ccccidx[i][3],ccccidx[i][4]);
        if (ccccidx[i][1]!=0&&ccccidx[i][2]&& \
            ccccidx[i][3]!=0&&ccccidx[i][4]!=0&&i<dim2) fprintf(out0, "+ \n");
    }
    fprintf(out0, ";");
    fclose(out0);
*/

    free(cc);
    free(cs);
    free(sc);
    free(ss);
    free(ccidx);

    return 0;
}

/* hypergeometric function 1F2 
   evaluated by truncating an infinite sum */
long double h1f2(long double a, long double b, long double c, long double z){
	int i, imax;
	long double err, tol, si, sii, f1f2;
	long double errabs, errrel;
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
long double LegendreQ(long double m, long double n, long double z){
    long double lq;
//    lq=pow(M_PI,.5)*exp(lgammal(m+n+.5))/(pow(2.,m+.5)*exp(lgammal(m+1.)));
    lq=pow(tanh(z),n)/pow(cosh(z),m+.5);
    lq=lq*h1f2(.5*(m+n+.5),.5*(m+n+1.5),m+1.,1./cosh(z)/cosh(z));
    return lq;
}

/* Derivative of Q_{m+1/2}^n (cosh(z)).
   Normalized to remove gamma function.
 */
long double DLegendreQ(int m, int n, long double z){
    long double dlq, lq1, lq2;
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
long double zetaf(long double rho, long double z, long double r0){
    long double zetax;
    zetax=atanhl(2.*rho*r0/(rho*rho+r0*r0+z*z));
    return zetax;
}

/* coordinate transformation r, z -> eta */
long double etaf(long double rho, long double z, long double r0){
    long double etax, xx;
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










