#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

#include <stdio.h>
#include "globheads.h"

#if defined(_SGI) || defined(_LINUX)
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#elif defined(_IBM)
#include <essl.h>
#define dnrm2   dnrm2
#define ddot    ddot
#define daxpy   daxpy
#define qsplit  qsplit
#define dscal   dscal
#define dgemv   dgemv
#define dgemm   dgemm
#define dgetrf  dgetrf
#define dgetri  dgetri
#define dgesvd  dgesvd
#define readmtc readmtc
#define csrcsc  csrcsc 
#define roscal  roscal
#define coscal  coscal
#define qsplit  qsplit
#else
#define dnrm2   dnrm2_
#define ddot    ddot_
#define daxpy   daxpy_
#define qsplit  qsplit_
#define dscal   dscal_
#define dgemv   dgemv_
#define dgemm   dgemm_
#define dgetrf  dgetrf_
#define dgetri  dgetri_
#define dgesvd  dgesvd_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#define qsplit  qsplit_
#endif

#define MAX_LINE        256
#define MAX_HBNAME      64
#define MAX_MAT			100



/* sets.c */


extern void *Malloc(int nbytes, char *msg); 
extern int setupCS(csptr amat, int len, int job); 
extern int cleanCS(csptr amat);
extern int nnz_cs (csptr A) ;
extern int cscpy(csptr amat, csptr bmat);

extern int setupVBMat(vbsptr vbmat, int n, int *nB);

extern int cleanVBMat(vbsptr vbmat); 
extern int nnzVBMat(vbsptr vbmat) ;
extern int memVBMat(vbsptr vbmat); 

extern void zrmC(int m, int n, BData data); 
extern void copyBData(int m, int n, BData dst, BData src, int isig); 
extern int CSRcs(int n, double *a, int *ja, int *ia, csptr mat, int 
		 rsa); 
extern int csrvbsrC(int job, int nBlk, int *nB, csptr csmat, vbsptr
		    vbmat);  
extern int col2vbcol( int col, vbsptr vbmat );
extern int COOcs(int n, int nnz,  double *a, int *ja, int *ia, csptr bmat);
void coocsr_(int*, int*, double*, int*, int*, double*, int*, int*);

/* MatOps.c */
extern int diag_scal(vbsptr vbmat);
extern int diagvec(vbsptr vbmat, BData x, BData y);
extern void matvec(csptr mata, double *x, double *y); 

extern void matvecz(csptr mata, double *x, double *y, double *z);
extern void vbmatvec(vbsptr vbmat, double *x, double *y);



extern int rpermC(csptr mat, int *perm); 
extern int cpermC(csptr mat, int *perm) ; 
extern int dpermC(csptr mat, int *perm) ; 
extern int CSparTran(csptr amat, csptr bmat, CompressType *compress);
extern double vbnorm2(int sz, double *a);
extern void Lsol(csptr mata, double *b, double *x);
extern void Usol(csptr mata, double *b, double *x);

















#ifndef _IBM
extern void dgemv(char*, int *, int*, double*, double *, int*,
		  double*, int*, double*, double*, int*);  
extern void dgemm(char*, char*, int*, int*, int*, double*, double*,
		  int*, double*, int*, double*, double*, int*) ;  
extern void dgetrf(int*, int*, double*, int*, int*, int*); 
extern void dgetri(int*, double*, int*, int*, double*,  int*, int* );
extern double dnrm2( int *, double *, int * );
extern void dscal(int*, double*, double*, int*); 
#endif 
extern int invGauss(int nn, double *A); 
extern int invSVD(int nn, double *A) ;

/* setblks.c */ 
extern int KeyComp(const void *vfst, const void *vsnd);
extern int init_blocks(csptr, int *, int **, int **, double, double*,double*);  

/* misc.c */
extern int SparTran(csptr amat, csptr bmat, int job, int flag); 
extern int coscalC(csptr mata, double *diag, int nrm);
extern void dscale(int n, double *dd, double *x, double * y);
extern void hilosort(csptr mat, int abval, int hilo);
extern void printmat(FILE *ft, csptr A, int i0, int i1);
extern void qqsort(int *ja, double *ma, int left, int right);
extern void qsort2C(int *ja, double *ma, int left, int right, int
		    abval); 
extern void qsort3i(int *wa, int *cor1, int *cor2, int left, int
		    right); 
extern void qsortC(int *ja, double *ma, int left, int right, int
		   abval); 
extern void qsortR2I(double *wa, int *cor1, int *cor2, int left, int
		     right); 
extern int qsplitC(double *a, int *ind, int n, int ncut);
extern int roscalC(csptr mata, double *diag, int nrm);
extern void swapj(int v[], int i, int j);
extern void swapm(double v[], int i, int j);

#endif 
