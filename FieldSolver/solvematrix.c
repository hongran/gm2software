/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
        #include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
        #include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include "lis.h"
LIS_INT main(LIS_INT argc, char* argv[])
{
    int mdim, ndim, id, jd, vid, vn, mmax, nmax;
    double tp;
    double tpl;
    FILE *mtrx, *vctr, *names, *coeffs, *vecdim;
    mtrx = fopen ("outmatrix.dat", "r");
    vecdim = fopen ("vecdim.txt", "r");
    fscanf (vecdim, "%d", &mdim);
    ndim=mdim;
    printf("%d x %d matrix \n", mdim, ndim);

    LIS_MATRIX A;
    LIS_VECTOR b,x;
    LIS_SOLVER solver;
    LIS_SCALAR tt, tv;
    LIS_INT my_rank;
    int int_nprocs,int_my_rank;
    LIS_INT err,i,j,n,gn,is,ie,iter;
    n  = ndim;
    lis_initialize(&argc, &argv);
//    tv = (LIS_SCALAR *)malloc( n*sizeof(LIS_SCALAR) );

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&int_nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&int_my_rank);
    my_rank = int_my_rank;
#else
    my_rank = 0;
#endif

    lis_matrix_create(LIS_COMM_WORLD,&A); 
//    lis_matrix_set_type(A, LIS_MATRIX_DNS); 
    err = lis_matrix_set_size(A,0,n);
    CHKERR(err);
    lis_matrix_get_size(A,&n,&gn);
    lis_matrix_get_range(A,&is,&ie);
    printf("matrix set (%d,%d)\n",is,ie);
    vn=0;
    for (i=0;i<mdim;i++) {
        for (j=0;j<mdim;j++) {
            //fscanf (mtrx, "%lg", &tpl);
	      fread(&tpl,sizeof(double),1,mtrx);
//            printf("%f\n",tp);
            if(i>=j) {
                tt=tpl;
                lis_matrix_set_value(LIS_INS_VALUE,i,j,tt,A);
                if(i!=j) lis_matrix_set_value(LIS_INS_VALUE,j,i,tt,A);
            }
        }
        vn++;
        if(vn==11) {
            printf("row %d loaded (%.2f % )\n",i,100.*((1.*i+1.)*(1.*i+1.))/(1.*ie)/(1.*ie));
            vn=1;
        }
    }
    printf("matrix loaded\n");
/*    for(i=is;i<ie;i++)
    {
        if( i>0   )  lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A);
        if( i<gn-1 ) lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A);
        lis_matrix_set_value(LIS_INS_VALUE,i,i,2.0,A);
    } */
//    lis_matrix_set_type(A,LIS_MATRIX_DNS);
    lis_matrix_set_type(A,LIS_MATRIX_CSR);
    lis_matrix_assemble(A);

    lis_vector_duplicate(A,&b);
    lis_vector_duplicate(A,&x);
    vctr = fopen ("outvector.dat", "r");
    for (i=0;i<mdim;i++) {
    //    fscanf (vctr, "%lg", &tpl);
        fread(&tpl,sizeof(double),1,vctr);
        tt=tpl;
        lis_vector_set_value(LIS_INS_VALUE,i,tt,b);
    }
    printf("vector loaded\n");
//    printf("start solving\n");
    lis_solver_create(&solver);
//    lis_solver_set_option("-print mem -i bicg -omp_num_threads 10 -maxiter 1000000",solver);
    lis_solver_set_option("-print out -i cg -omp_num_threads 10 -maxiter 1000000 -tol 1.0e-9",solver);
    lis_solver_set_optionC(solver);
    lis_solve(A,b,x,solver);
    lis_solver_get_iter(solver,&iter);
    if (my_rank==0)
      {
#ifdef _LONG__LONG
	printf("number of iterations = %lld\n",iter);
#else
	printf("number of iterations = %d\n",iter);
#endif
	printf("\n");
      }
//    lis_vector_print(x);
    printf("writing to disk\n");
    names = fopen ("outnames.txt", "r");
    coeffs = fopen ("fitcoeffN.txt", "w");
    fscanf (names, "%d", &mmax);
    fscanf (names, "%d", &nmax);
    fprintf (coeffs, "%d\n", mmax);
    fprintf (coeffs, "%d\n", nmax);
    printf("writing to disk\n");
    lis_output_vector(x, LIS_FMT_PLAIN, "fitc.txt");
//    lis_vector_get_values(x,0,ndim,tv);
/*
    for (j=0;j<ndim;j++) {
        fscanf (names, "%d \t %d \t %d", &vn, &id, &jd);
        lis_vector_get_value(x,j,tv);
        printf("ee\n");
        tpl=tv[j];
        printf("%d, %Lf \n", j+1,tpl);
        fprintf (coeffs, "%d \t %d \t %d \t %.20Le\n", \
                          vn, id, jd, tpl);
    } */
                                                                                    printf("wrote to disk\n");


    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_solver_destroy(solver);
    lis_finalize();
    return 0;
}
