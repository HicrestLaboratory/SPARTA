#include "mkl_utilities.h"
#include "globheads.h"
#include "protos.h"
#include "mkl.h"
#include <stdio.h>
#include <iostream>

using namespace std;

//polymorphism for mkl_gemm_custom;
//----------------------------------------------------------------------------------
//FLOAT
void mkl_gemm_custom(const CBLAS_LAYOUT Layout,
                      const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                      const MKL_INT m, const MKL_INT n, const MKL_INT k,
                      const float alpha, const float *a, const MKL_INT lda,
                      const float *b, const MKL_INT ldb,
                     const float beta, float *c, const MKL_INT ldc){
    
    cblas_sgemm (Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

    }

//DOUBLE
void mkl_gemm_custom(const CBLAS_LAYOUT Layout,
                  const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                  const MKL_INT m, const MKL_INT n, const MKL_INT k,
                  const double alpha, const double *a, const MKL_INT lda,
                  const double *b, const MKL_INT ldb,
                     const double beta, double *c, const MKL_INT ldc){

    cblas_dgemm (Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

    }

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//polymorphism for mkl_spmm_custom;
//----------------------------------------------------------------------------------
//FLOAT
void mkl_spmm_custom(const sparse_operation_t operation,
                     const float alpha, const sparse_matrix_t A,
                     const struct matrix_descr descr,
                     const sparse_layout_t layout,
                     const float *x, const MKL_INT columns,
                     const MKL_INT ldx,
                     const float beta,
                     float *y, const MKL_INT ldy){
    
    mkl_sparse_s_mm (operation,
                     alpha, A,
                     descr,
                     layout,
                     x, columns, ldx,
                     beta,
                     y, ldy);

    }

//DOUBLE
void mkl_spmm_custom(const sparse_operation_t operation,
                 const double alpha, const sparse_matrix_t A,
                 const struct matrix_descr descr,
                 const sparse_layout_t layout,
                 const double *x, const MKL_INT columns,
                 const MKL_INT ldx,
                 const double beta,
                 double *y, const MKL_INT ldy){

    mkl_sparse_d_mm (operation,
                 alpha, A,
                 descr,
                 layout,
                 x, columns, ldx,
                 beta,
                 y, ldy);

    }
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------





//polymorphism for mkl_batch_custom;
//----------------------------------------------------------------------------------
void mkl_batch_custom(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE* transa_array, const CBLAS_TRANSPOSE* transb_array,
            const MKL_INT* m_array, const MKL_INT* n_array, const MKL_INT* k_array,
            const float* alpha_array,
            const float **a_array, const MKL_INT* lda_array,
            const float **b_array, const MKL_INT* ldb_array,
            const float* beta_array,
            float **c_array, const MKL_INT* ldc_array,
            const MKL_INT group_count, const MKL_INT* group_size){
    
    cblas_sgemm_batch (Layout, transa_array, transb_array, m_array,  n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array,  beta_array, c_array, ldc_array, group_count, group_size);

    }

void mkl_batch_custom(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE* transa_array, const CBLAS_TRANSPOSE* transb_array,
            const MKL_INT* m_array, const MKL_INT* n_array, const MKL_INT* k_array,
            const double* alpha_array,
            const double **a_array, const MKL_INT* lda_array,
            const double **b_array, const MKL_INT* ldb_array,
            const double* beta_array,
            double **c_array, const MKL_INT* ldc_array,
            const MKL_INT group_count, const MKL_INT* group_size){
    
    cblas_dgemm_batch (Layout, transa_array, transb_array, m_array,  n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array,  beta_array, c_array, ldc_array, group_count, group_size);

    }
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


//polymorphism for mkl_create_csr_custom;
//----------------------------------------------------------------------------------
//FLOAT
void mkl_create_csr_custom(sparse_matrix_t *A,
                           const sparse_index_base_t indexing,
                           const MKL_INT rows, const MKL_INT cols,
                           MKL_INT *rows_start, MKL_INT *rows_end,
                           MKL_INT *col_indx,
                           float *values){
    
    mkl_sparse_s_create_csr(A, indexing, rows, cols, rows_start, rows_end, col_indx, values);
    }

//DOUBLE
void mkl_create_csr_custom(sparse_matrix_t *A,
                           const sparse_index_base_t indexing,
                           const MKL_INT rows, const MKL_INT cols,
                           MKL_INT *rows_start, MKL_INT *rows_end,
                           MKL_INT *col_indx,
                           double *values){

    mkl_sparse_d_create_csr(A, indexing, rows, cols, rows_start, rows_end, col_indx, values);
}


void convert_to_MKL(SparMat &spmt, sparse_matrix_t &A){
    
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
    
    int n = spmt.n;
    MKL_INT cols = n;
    MKL_INT rows = n;
    MKL_INT* rows_start= new MKL_INT[n];
    MKL_INT* rows_end = new MKL_INT[n];
    
    rows_start[0] = 0;
    rows_end[0] = spmt.nzcount[0];
    for (int i = 0; i<n;i++){
        rows_start[i] = spmt.nzcount[i-1] + rows_start[i - 1];
        rows_end[i] = spmt.nzcount[i] + rows_end[i - 1];
    }

    
    int nzs = rows_end[n-1];
    DataT* values = new DataT[nzs];
    MKL_INT* col_indx = new MKL_INT[nzs];
    
    int j = 0;
    for (int row = 0; row < n; row++){
        for (int i = 0; i <spmt.nzcount[row]; i++){
            col_indx[j] = spmt.ja[row][i];
            j++;
        }
    }
    
    j = 0;
    for (int row = 0; row < n; row++){
        for (int i = 0; i <spmt.nzcount[row]; i++){
            values[j] = spmt.ma[row][i];
            j++;
        }
    }
    
    mkl_create_csr_custom (&A, indexing, rows, cols, rows_start,  rows_end, col_indx, values);

}



//BATCH multiply a n-by-n block matrix VBMat by a (column major) n-by-k matrix X.
//store result in (already initialized) (column major) n-by-k matrix Y;
void mkl_blockmat_batch_multiply(const VBSparMat &VBMat, DataT *X, int X_cols, DataT *Y){
    int N = VBMat.n, *bsz = VBMat.bsz;
    int Lsz,Hsz,col;
    int mat_n = bsz[N];


    int h_scan = 0; //horizontal scan counter
    int batch_count  = -1; //counter for gemm in a single batch
    
    MKL_INT     ms[N];
    MKL_INT    ns[N];
    MKL_INT    ks[N];
    MKL_INT    lda_array[N];
    MKL_INT     ldb_array[N];
    MKL_INT    ldc_array[N];

    CBLAS_TRANSPOSE    transA[N];
    CBLAS_TRANSPOSE    transB[N];

    DataT    alpha[N];
    DataT    beta[N];

    DataT *a_array[N];
    DataT *b_array[N];
    DataT *c_array[N];

    MKL_INT    size_per_grp[N];
    

    while (h_scan < N  & batch_count != 0){ //exit when there are no more block_columns to process
            
        //loop vertically through block rows
        batch_count = 0;
        for(int i = 0; i < N; i++ ) {
            Hsz = bsz[i+1] - bsz[i];

            //loop horizontaly through block columns
            if (h_scan < VBMat.nzcount[i]){

                batch_count++;
                
                col = VBMat.ja[i][h_scan];
                Lsz = bsz[col+1] - bsz[col];
                //multiply the block by the matrix
                //define the sub-matrices
                DataT* block = (VBMat.ba)[i][h_scan]; //access block i,j in column major order.
                DataT* blockY = Y + bsz[i];     //i indicates the vertical block of Y that is going to be updated
                DataT* blockX = X + bsz[col];     //col indicates the vertical block of X that is going to be multiplied with the (i,j)block of VBMat
                
                ms[batch_count-1] = Hsz;
                ns[batch_count-1] = X_cols;
                ks[batch_count-1] = Lsz;
                
                lda_array[batch_count-1] = Hsz;
                ldb_array[batch_count-1] = mat_n;
                ldc_array[batch_count-1] = mat_n;

                transA[batch_count-1] = CblasNoTrans;
                transB[batch_count-1] = CblasNoTrans;

                alpha[batch_count-1] = 1.;
                beta[batch_count-1] = 1.;

                a_array[batch_count-1] = block;
                b_array[batch_count-1] = blockX;
                c_array[batch_count-1] = blockY;

                size_per_grp[batch_count-1] = 1;
                
            }
        
        }
	h_scan++;
        if(batch_count > 0) {
            mkl_batch_custom (CblasColMajor, transA, transB,
                              ms, ns, ks, alpha,
                              (const DataT **) a_array, lda_array,
                              (const DataT**) b_array, ldb_array,
                              beta,
                              c_array, ldc_array,
                              batch_count, size_per_grp);
        }
    }
 
}

//TODO fix bug that alters some rows from the right result
//multiply a n-by-n block matrix VBMat by a (column major) n-by-k matrix X.
//store result in (already initialized) (column major) n-by-k matrix Y;
void mkl_blockmat_multiply(const VBSparMat &VBMat, DataT *X, int X_cols, DataT *Y){
    int N = VBMat.n, *bsz = VBMat.bsz;
    int Lsz,Hsz,col;
    int mat_n = bsz[N];
    DataT alpha = 1.0;
    DataT beta = 1.0;
   
    //loop vertically through block rows
    for(int i = 0; i < N; i++ ) {
        Hsz = bsz[i+1] - bsz[i];

        //loop horizontaly through block columns
        for(int j = 0; j<VBMat.nzcount[i]; j++){
            col = VBMat.ja[i][j];
            Lsz = bsz[col+1] - bsz[col];
            //multiply the block by the matrix
            //define the sub-matrices
            const DataT* block = (VBMat.ba)[i][j];  //access block i,j in column major order.
            DataT* blockY = Y + bsz[i];             //i indicates the vertical block of Y that is going to be updated
            const DataT* blockX = X + bsz[col];     //col indicates the vertical block of X that is going to be multiplied with the (i,j)block of VBMat
            mkl_gemm_custom(CblasColMajor, CblasNoTrans, CblasNoTrans, Hsz, X_cols, Lsz, alpha, block, Hsz, blockX, mat_n, beta, blockY, mat_n);
            }
        }
}

