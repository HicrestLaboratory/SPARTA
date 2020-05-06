#ifndef __VBLOCK_HEADER_H__
#define __VBLOCK_HEADER_H__

#define MAX_BLOCK_SIZE   100

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define B_DIM(bs,i)      (bs[i+1]-bs[i])

typedef float DataT; //precision for matrix entries
typedef DataT *BData;

typedef struct SpaFmt {
/*--------------------------------------------- 
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
  int n;
  int *nzcount;  /* length of each row */
  int **ja;      /* pointer-to-pointer to store column indices  */
  DataT **ma;   /* pointer-to-pointer to store nonzero entries */
    
  
} SparMat, *csptr;


typedef struct VBSpaFmt {
    int n;        /* the block row dimension of the matrix      */
    int *bsz;     /* the row/col of the first element of each   */
                  /* diagonal block                             */
    int *nzcount;  /* length of each row                         */
    int **ja;     /* pointer-to-pointer to store column indices */
    BData **ba;   /* pointer-to-pointer to store nonzero blocks */
    BData *D;     /* to store inversion of diagonals            */
} VBSparMat, *vbsptr;

typedef struct __CompressType
{
  int grp;   /* -1: begin new group, >=0: belong to grp-th row */
  int count; /* block size, valid only if grp = -1 */
} CompressType;


#endif  /* __VBLOCK_HEADER_H__ */

