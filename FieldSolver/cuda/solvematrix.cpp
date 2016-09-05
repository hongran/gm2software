/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
 * This sample implements a conjugate gradient solver on GPU
 * using CUBLAS and CUSPARSE
 *
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

// Utilities and system includes
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <helper_cuda.h>       // helper function CUDA error checking and initialization

const char *sSDKname     = "conjugateGradient";


int main(int argc, char **argv)
{
    int M = 0, N = 0, nz = 0, *I = NULL, *J = NULL;
    double *val = NULL;
    const double tol = 1e-9;
    const int max_iter = 1000000;
    double *x;
    double *rhs;
    double a, b, na, r0, r1;
    int *d_col, *d_row;
    double *d_val, *d_x, dot;
    double *d_r, *d_p, *d_Ax;
    int k;
    double alpha, beta, alpham1;
    //file pointers
    FILE *mtrx, *vctr, *names, *coeffs, *vecdim , *outvec;

    // This will pick the best possible CUDA capable device
    cudaDeviceProp deviceProp;
    int devID = findCudaDevice(argc, (const char **)argv);

    if (devID < 0)
    {
        printf("exiting...\n");
        exit(EXIT_SUCCESS);
    }

    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

    // Statistics about the GPU device
    printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
           deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);

    int version = (deviceProp.major * 0x10 + deviceProp.minor);

    if (version < 0x11)
    {
        printf("%s: requires a minimum CUDA compute 1.1 capability\n", sSDKname);

        // cudaDeviceReset causes the driver to clean up all state. While
        // not mandatory in normal operation, it is good practice.  It is also
        // needed to ensure correct operation when the application is being
        // profiled. Calling cudaDeviceReset causes all profile data to be
        // flushed before the application exits
        cudaDeviceReset();
        exit(EXIT_SUCCESS);
    }

    ///From original code///////////////////////////////
    vecdim = fopen ("vecdim.txt", "r");
    fscanf (vecdim, "%d", &N);
    fclose(vecdim);
    M=N;
    nz=M*N;
    printf("%d x %d matrix \n", M, N);

    I = (int *)malloc(sizeof(int)*(N+1));
    J = (int *)malloc(sizeof(int)*nz);
    val = (double *)malloc(sizeof(double)*nz);
    /* load matrix in CSR format */
    //Loading matrix
    memset(val, 0, sizeof(double)*M*N);
    //mtrx = fopen ("outmatrix.txt", "r");
    mtrx = fopen ("outmatrix.dat", "r");
    int vn=0;
    double tpl;
    int Jdx = 0;
    for (int i=0;i<M;i++) {
      for (int j=0;j<N;j++) {
	//fscanf (mtrx, "%lg", &tpl);
	fread(&tpl,sizeof(double),1,mtrx);
	//            printf("%f\n",tp);
	if(i>=j) {
	  val[i*N+j]=tpl;
	  if(i!=j) val[j*M+i]=tpl;
	}
	J[Jdx]=j;
	Jdx++;
      }
      I[i] = i*N;
      vn++;
      if(vn==11) {
	printf("row %d loaded \n",i);
	vn=1;
      }
    }
    I[M]=M*N;
    fclose(mtrx);
    printf("matrix loaded\n");

    /////////////////////////////////////////////////////

    /* Generate a random tridiagonal symmetric matrix in CSR format */
/*    M = N = 1048576;
    nz = (N-2)*3 + 4;
    I = (int *)malloc(sizeof(int)*(N+1));
    J = (int *)malloc(sizeof(int)*nz);
    val = (double *)malloc(sizeof(double)*nz);
    genTridiag(I, J, val, N, nz);
*/
    x = (double *)malloc(sizeof(double)*N);
    rhs = (double *)malloc(sizeof(double)*N);

    for (int i = 0; i < N; i++)
    {
        rhs[i] = 0.0;
        x[i] = 0.0;
    }
    //Load vector
    //vctr = fopen ("outvector.txt", "r");
    vctr = fopen ("outvector.dat", "r");
    for (int i=0;i<N;i++) {
      //fscanf (vctr, "%lg", &tpl);
      fread(&tpl,sizeof(double),1,vctr);
      rhs[i]=tpl;
    }
    fclose(vctr);
    printf("vector loaded\n");


    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&cublasHandle);

    checkCudaErrors(cublasStatus);

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);

    checkCudaErrors(cusparseStatus);

    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);

    checkCudaErrors(cusparseStatus);

    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    checkCudaErrors(cudaMalloc((void **)&d_col, nz*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_val, nz*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_x, N*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_r, N*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(double)));
    checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(double)));

    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, rhs, N*sizeof(double), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasDscal(cublasHandle, N, &b, d_p, 1);
            cublasStatus = cublasDaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublasStatus = cublasDdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasDaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasDaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaThreadSynchronize();
        printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }

    cudaMemcpy(x, d_x, N*sizeof(double), cudaMemcpyDeviceToHost);

    double rsum, diff, err = 0.0;

    for (int i = 0; i < N; i++)
    {
        rsum = 0.0;

        for (int j = I[i]; j < I[i+1]; j++)
        {
            rsum += val[j]*x[J[j]];
        }

        diff = fabs(rsum - rhs[i]);

        if (diff > err)
        {
            err = diff;
        }
    }

    //Output
    names = fopen ("outnames.txt", "r");
    coeffs = fopen ("fitcoeffN.txt", "w");
    int mmax,nmax;
    fscanf (names, "%d", &mmax);
    fscanf (names, "%d", &nmax);
    fprintf (coeffs, "%d\n", mmax);
    fprintf (coeffs, "%d\n", nmax);
    fclose(names);
    fclose(coeffs);
    printf("writing to disk\n");
    outvec = fopen("fitc.txt","w");
    for (int i=0;i<M;i++) {
      fprintf(outvec,"%.17g\n",x[i]);
    }
    fclose(outvec);


    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    free(I);
    free(J);
    free(val);
    free(x);
    free(rhs);
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ax);

    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();

    printf("Test Summary:  Error amount = %f\n", err);
    exit((k <= max_iter) ? 0 : 1);
}
