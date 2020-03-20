
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "globheads.h"
#include "protos.h"
using namespace std;

typedef struct __KeyType
{
  int var;   /* row number */
  int key;   /* hash value */
} KeyType;

int KeyComp( const void *vfst, const void *vsnd )
{
  KeyType *fst = (KeyType *)vfst, *snd = (KeyType *)vsnd;
  if( fst->key == snd->key ) {
    if( fst->var < snd->var )
      return -1;
    return 1;
  }
  if( fst->key < snd->key )
    return -1;
  return 1;
}

int setupCS(csptr amat, int len, int job)
{

/*----------------------------------------------------------------------
| Initialize SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|     len   =  size of matrix
|     job   =  0: pattern only
|              1: data and pattern
|
| On return:
|===========
|
|  amat->n
|      ->*nzcount
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
    
   amat->n = len;
   amat->nzcount = new int[len];

   amat->ja = (int **) malloc( len*sizeof(int *));


   if( job == 1 )
       amat->ma = (DataT **) malloc( len*sizeof(DataT *));
   else
       amat->ma = NULL;

   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCS
|--------------------------------------------------------------------*/

int cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|--------------------------------------------------------------------*/
   /*   */
  int i;
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
  for (i=0; i<amat->n; i++) {
    if (amat->nzcount[i] > 0) {
      if( amat->ma ) free(amat->ma[i]);
      free(amat->ja[i]);
    }
  }
  if (amat->ma) free(amat->ma);
  free(amat->ja);
  free(amat->nzcount);
  free(amat);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/


int CSparTran( csptr amat, csptr bmat, CompressType *compress )
{
/*----------------------------------------------------------------------
| Finds the compressed transpose of a matrix stored in SpaFmt format.
| Patterns only.
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat)     = a matrix stored in SpaFmt format.
| (compress) = quotient graph of matrix amat
|
| on return:
| ----------
| (bmat)     = the compressed transpose of (mata) stored in SpaFmt
|              format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/

   int i, j, *ind, nzcount, pos, size=amat->n, *aja;
   ind = bmat->nzcount;


    for (i=0; i<size; i++)
      ind[i] = 0;
/*-------------------- compute lengths  */
   for (i=0; i<size; i++) {
     if( compress[i].grp != -1 ) continue;
     aja = amat->ja[i];

     nzcount = amat->nzcount[i];

     for (j=0; j < nzcount; j++) {
       pos = aja[j];

       if( compress[pos].grp == -1 ) {
         ind[pos]++;
       }
     }
   }

/*--------------------  allocate space  */
   for (i=0; i<size; i++) {
      if( ind[i] == 0 ) {
        bmat->ja[i] = NULL;
        continue;
      }
      bmat->ja[i] = (int *)malloc( ind[i]*sizeof(int));
      ind[i] = 0; /* indicate next available position of each row */
   }
/*--------------------  now do the actual copying  */
   for (i=0; i<size; i++) {
      if( compress[i].grp != -1 ) continue;
      aja = amat->ja[i];
      nzcount = amat->nzcount[i];
      for (j = 0; j < nzcount; j++) {
         pos = aja[j];
         if( compress[pos].grp == -1 ) {
           bmat->ja[pos][ind[pos]] = i;
           ind[pos]++;
         }
      }
   }
   return 0;
}

int csrvbsrC( int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat )
{
/*----------------------------------------------------------------------
 *  Compressed C-style Sparse Row to C-style Various Block Sparse Row
 *----------------------------------------------------------------------
 *
 * This  subroutine converts a matrix stored  in a C-style SpaFmt format
 * into a C-style various block SpaFmt format
 *
 * NOTE: the initial matrix does not have to have a block structure.
 * zero padding is done for general sparse matrices.
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 * job   = if job == 0 on entry, pattern only is generated.
 *
 * nBlk  = integer equal to the dimension of block matrix.
 *
 * nB    = integer array of diagonals' block size
 *
 * csmat = Sparse Row format Matrix
 *
 * on return:
 *-----------
 *
 * vbmat = Various Block Sparse Row format Matrix
 *
 * ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int n, i, j, k;
    int nnz, szofBlock, ipos, b_row, b_col, br, bc;
    int *iw = NULL;
    
    int counter = 0;
        
    n  = csmat->n;            /* size of the original matrix          */
    setupVBMat( vbmat, nBlk, nB );
    iw = (int *)malloc(sizeof(int)*nBlk);
    for( i = 0; i < nBlk; i++ ) iw[i] = 0;
    b_row = -1;
    for( i = 0; i < n; i += nB[b_row] ) {
        vbmat->nzcount[++b_row] = 0;
        
        /* calculate nzcount of the (b_row)-th row of the block matrix */
        for( j = i; j < i+nB[b_row]; j++ ) {
            
            int nnz_j = csmat->nzcount[j];
            for( k = 0; k < nnz_j; k++ ) {
                /* get the column ID of block matrix by giving the column ID
                   of the original matrix */
                b_col = col2vbcol( csmat->ja[j][k], vbmat );
                if( iw[b_col] == 0 ) {
                    iw[b_col] = 1;
                    vbmat->nzcount[b_row]++;
                }
            }
        }
        if( 0 == ( nnz = vbmat->nzcount[b_row] ) ) continue;
        vbmat->ja[b_row] = (int *)malloc( sizeof(int)*nnz);

        /* calculate the pattern of the (b_row)-th row of the block matrix */
        for( j = 0, ipos = 0; j < nBlk; j++ ) {
            if( iw[j] != 0 ) {
                vbmat->ja[b_row][ipos] = j;
                iw[j] = ipos;
                ipos++;
            }
        }
            
        
        if( job == 0 ) goto NEXT_ROW;  /* stop here if patterns only */

        /* copy data to the (b_row)-th row of the block matrix from the
           original matrix */
        vbmat->ba[b_row] = (BData *)malloc( sizeof(BData)*nnz);
        for( j = 0; j < nnz; j++ ) {
            szofBlock = sizeof(DataT)*nB[b_row]*nB[vbmat->ja[b_row][j]];
            vbmat->ba[b_row][j] = (BData)malloc( szofBlock);
        }
        
        for( j = i; j < i+nB[b_row]; j++ ) {
            for( k = 0; k < csmat->nzcount[j]; k++ ) {
                /* get the column ID of block matrix by giving the column ID
                   of the original matrix */
                b_col = col2vbcol( csmat->ja[j][k], vbmat );
                ipos = iw[b_col];
                br = j - i;
                bc = csmat->ja[j][k] - vbmat->bsz[b_col];
                DATA(vbmat->ba[b_row][ipos],nB[b_row],br,bc) = csmat->ma[j][k];
            }
        }
NEXT_ROW:
        /* reset iw */
        for( j = 0; j < nnz; j++ ) iw[vbmat->ja[b_row][j]] = 0;
    }
    
    free( iw );
    return 0;
}

int col2vbcol( int col, vbsptr vbmat )
{
/*---------------------------------------------------------------------
 * get the column ID of block matrix by giving the column ID of the original
 * matrix
 *--------------------------------------------------------------------*/
    int *bsz = vbmat->bsz, n = vbmat->n;
    int begin = 0, mid, end = n-1;
    while( end - begin > 1 ) {
        mid = (begin+end)/2;
        if( col < bsz[mid] ) {
            end = mid;
        } else if( col >= bsz[mid+1] ) {
            begin = mid;
        } else {
            return mid;
        }
    }
    if( col >= bsz[end] ) {
        return end;
    }
    return begin;
}

int setupVBMat( vbsptr vbmat, int n, int *nB )
{
/*----------------------------------------------------------------------
| Initialize VBSpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat ) =  Pointer to a VBSpaFmt struct.
|       n   =  size of block matrix
|      nB   =  size of diagonal block, so the real size of the matrix
|              is nB[0] + nB[1] + ... + nB[n-1]
|              do nothing if nB is NULL
|
| On return:
|===========
|
| vbmat->n
|      ->*bsz
|      ->*nzcount
|      ->**ja
|      ->**ba
|      ->*D
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int i;
    vbmat->n = n;
    if( nB ) {
        vbmat->bsz = (int *)malloc( sizeof(int)*(n+1));
        vbmat->bsz[0] = 0;
        for( i = 1; i <= n; i++ ) {
            vbmat->bsz[i] = vbmat->bsz[i-1] + nB[i-1];
        }
    } else
        vbmat->bsz = NULL;
    vbmat->nzcount = (int *)malloc( sizeof(int)*n);
    vbmat->ja = (int **)malloc( sizeof(int *)*n);
    vbmat->ba = (BData **)malloc( sizeof(BData *) * n);
    vbmat->D = NULL;
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupVBMat
|--------------------------------------------------------------------*/

int cleanVBMat( vbsptr vbmat )
{
/*----------------------------------------------------------------------
| Free up memory allocated for VBSpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a VBSpaFmt struct.
|--------------------------------------------------------------------*/
    int i, j;
    if( vbmat == NULL ) return 0;
    if( vbmat->n < 1 ) return 0;
    
    for( i = 0; i < vbmat->n; i++ ) {
        if( vbmat->nzcount[i] > 0 ) {
            free( vbmat->ja[i] );
            if( vbmat->ba && vbmat->ba[i] ) {
                for( j = 0; j < vbmat->nzcount[i]; j++ ) {
                    free( vbmat->ba[i][j] );
                }
                free( vbmat->ba[i] );
            }
        }
        if( vbmat->D && vbmat->D[i] ) free( vbmat->D[i] );
    }
    if( vbmat->D ) free( vbmat->D );
    free( vbmat->ja );
    if( vbmat->ba ) free( vbmat->ba );
    free( vbmat->nzcount );
    if( vbmat->bsz ) free( vbmat->bsz );
    free( vbmat );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBMat
|--------------------------------------------------------------------*/



int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm,
                 double eps, double *t_hash, double *t_angle )
{
/*----------------------------------------------------------------------------
 * Setup Blocks ( rows and columns might be permuted to get better results )
 *----------------------------------------------------------------------------
 * Na Li, Aug 2001
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * csmat   = a matrix stored in SpaFmt format
 * eps     = parameter for deciding when to do a union of two rows
 *           into the same group.  Two rows u and v are merged into a
 *           block  when cos(<u,v>) == (u,v)/(|u|*|v|), is > eps.
 *           eps should be <= 1.
 *----------------------------------------------------------------------------
 * on return:
 * ==========
 * csmat   = matrix stored in SpaFmt format after permutation
 * pnBlock = dimension of the block matrix
 * pnB     = dimension of each block
 *
 *----------------------------------------------------------------------------
 * Combination of hash method and angle method:
 *----------------------------------------------------------------------------
 * Designed for the matrices with symmetric patterns
 * (1) Hash method
 *     a. Calculate hash values
 *     b. qsort rows according to their hash values
 *     c. Get compressed graph as the following format:
 * (2) Angle method
 *     a. Calculate A^T
 *     b. for i-th row, calculate dot product (row_i, row_j) using A*A^T
 *        algorithm where j = i+1, ..., n-1 and group[j] == -1
 *        if cos( <row_i, row_j> ) = (row_i,row_j)/|row_i||row_j| is > eps,
 *        we merge row_i and row_j by resetting
 *        group[j] = i and size[i] = size[i]+size[j]
 *--------------------------------------------------------------------------*/
  int n = csmat->n, nBlock = 0, i, j, k;
  csptr at = NULL;
  KeyType *group = NULL;
  CompressType *compress = NULL;
  int *perm = NULL, *nB = NULL;
  int nzcount0, nzcount, key0, key, *ja0, *ja, row0, row, newblock;
  int *iw = NULL, *jbuf = NULL;
  int cnt, pos, nnz_i, row_j, col, bkcnt;
  int nextBlockID, nextBlockPos, belongTo, grp;
  double eps_2 = eps * eps, t1, t2;

  group = (KeyType *)malloc( n*sizeof(KeyType));
  compress = (CompressType *)malloc( n*sizeof(CompressType));
  perm = (int *)malloc( n * sizeof(int));
  iw = perm; /* iw and perm array can share memory here because they will
          * never be used at the same time */
  for( i = 0; i < n; i++ ) {
    iw[i] = 0;
    compress[i].grp = -1;
  }
/*-------------------- compress matrix based on hash algorithm */
/*-------------------- get hash value of each row */
  for( i = 0; i < n; i++ ) {
    nzcount = csmat->nzcount[i];
    key = 0;
    ja = csmat->ja[i];
    for( j = 0; j < nzcount; j++ )
      key += ja[j]+1;
    group[i].key = key;
    group[i].var = i;
  }
/*-------------------- sort rows -- uses function KeyComp */
  qsort( group, n, sizeof(KeyType), KeyComp );

/*-------------------- compress matrix */
  for( i = 0; i < n; i++ ) {
    row0 = group[i].var;
    if( compress[row0].grp != -1 ) continue; /* already assigned */
    key0 = group[i].key;
    nzcount0 = csmat->nzcount[row0];
    ja0 = csmat->ja[row0];
/*-------------------- beginning of new block. set .grp and .count */
    compress[row0].grp = -1;
    compress[row0].count = 1;
/*-------------------- loop over all rows having same check-sum keys */
    for( j = i + 1; j < n; j++ ) {
      key = group[j].key;
      if( key != key0 ) break;
      row = group[j].var;
      if( compress[row].grp != -1 ) continue; /* already assigned */
      nzcount = csmat->nzcount[row];
      if( nzcount != nzcount0 ) continue;
      ja = csmat->ja[row];
      newblock = 0;
/*-------------------- compare patterns of the rows             */
      for( k = 0; k < nzcount; k++ ) iw[ ja0[k] ] = 1;
      for( k = 0; k < nzcount; k++ ) {
        if( iw[ ja[k] ] == 0 ) {
          newblock = 1;
          break;
        }
      }
      for( k = 0; k < nzcount; k++ ) iw[ ja0[k] ] = 0; /* reset iw */
/*-------------------- row belongs to group row0                    */
      if( !newblock ) {
        compress[row].grp = row0;
        compress[row0].count++;
      }
    }
  }
    
  nB = (int *)malloc( n * sizeof(int));
  jbuf = (int *)malloc( n * sizeof(int));

/*-------------------- compress matrix based on angle algorithm */
/*-------------------- calculate compressed A^T                 */
  at = (csptr)malloc( sizeof(SparMat));
  setupCS( at, n, 0 );
  if( CSparTran( csmat, at, compress ) != 0 )
    return -1;

/*----------------------------------------------------------------------------
 * only the row representing beginning of block satisfies:
 *    compress[row].grp = -1, so far.
 * how many such rows is up to the compression rate of hash compression
 * algorithm we did above
 *--------------------------------------------------------------------------*/

/*---------------------------------------------------------------
 * use group to backup original compressed matrix by Hash method.
 * It is very important because the array 'compress' will be changed
 * during Angle method. Or, we'll get incorrect inner product.
 *--------------------------------------------------------------*/
  for( i = 0; i < n; i++ ) {
    group[i].var = compress[i].grp;
    group[i].key = compress[i].count;
  }

  for( i = 0; i < n; i++ ) {
    if( compress[i].grp != -1 ) continue;
    nB[nBlock] = compress[i].count; /* !!! not 1 here */
    cnt = 0;
/*-------------------- calculate (u,v_j ), j = i+1,...,n, using product
 *-------------------- algorithm of A * A_T */
    nnz_i = csmat->nzcount[i];
    for( j = 0; j < nnz_i; j++ ) {
      row_j = csmat->ja[i][j];
      if( group[row_j].var != -1 ) /* i.e. original compress[row_j].grp */
        continue;
      bkcnt = group[row_j].key;    /* i.e. original compress[row_j].count */
      for( k = at->nzcount[row_j] - 1; k >= 0; k-- ) {
    col = at->ja[row_j][k];
    if( col <= i ) break;
        if( compress[col].grp != -1 ) continue; /* needed because compress
                                           array is dynamically updated */
    if( iw[col] == 0 ) { /* new nonzero of (u,v_j) */
      jbuf[cnt] = col;
      cnt++;
    }
        iw[col] += bkcnt; /* correct for matrix with symmetric pattern */
      }
    }
/*-------------------- set group for row i and reset iw */
    for( j = 0; j < cnt; j++ ) {
      pos = jbuf[j];
      if( iw[pos] * iw[pos] >= eps_2 * nnz_i * csmat->nzcount[pos] ) {
		compress[pos].grp = i;
		nB[nBlock] += compress[pos].count; /* !!! not 1 here */
      }
      iw[pos] = 0; /* reset iw */
    }
    nBlock++; /* begin new block, add block count by 1 */
  } /* end loop i */

/*-------------------- free group                                   */
  if( group ) {
    /* no need group array any more */
    free( group );
    group = NULL;
  }

  *pnBlock = nBlock;
  *pnB = (int *)malloc( nBlock * sizeof(int));
  for( i = 0; i < nBlock; i++ ) {
    if( nB[i] > MAX_BLOCK_SIZE ) {
      fprintf( stderr, "Block of size = %d exceeds MAX_BLOCK_SIZE\n", nB[i] );
      return -1;
    }
    (*pnB)[i] = nB[i];
  }

/*-------------------- calculate permutation array -  Array nB will
 * be used to store next  available position in each  block */
  nextBlockID = 0;
  nextBlockPos = 0;
  for( i = 0; i < n; i++ ) {
    if( compress[i].grp == -1 ) {
      perm[i] = nextBlockPos;
      nextBlockPos += (*pnB)[nextBlockID++];
      nB[i] = 1;
    } else {
      belongTo = compress[i].grp;
      grp = compress[belongTo].grp;
      if( grp != -1 ) /* should find the final beginning of block */
    belongTo = grp;
      perm[i] = perm[belongTo] + nB[belongTo];
      nB[belongTo]++;
    }
  }
  *pperm = perm;

  cleanCS( at );
  free( nB );
  free( jbuf );
  free( compress );

  return 0;
}
