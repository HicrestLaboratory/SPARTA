# SPARTA
**SPAR**se acceleration on **T**ensor **A**rchitecture

The project aims to investigate new data structures and compression algorithms for exploiting new architecture capabilities, specifically designed for deep learning, to accelerate **sparse and irregular** applications, such as graph analytics and arbitrary sparse DNN and GCN. SPARTA also looks at productivity and performance portability across different AI accelerators by providing an abstraction layer.  

The repository contains stable code for reordering and compressing sparse matrices into dense block data-structures.
The reordering algorithm matches rows (or columns) with similar patterns and builds dense blocks. 
The similarity of patterns is first determined with a hash function, and then refined with a tunable algorithm, which matches patterns with high cosine similarity.

The repository also contains code for sparse-dense matrix-matrix multiplication that exploits the dense block data-structure.

Input sparse matrices are stored in Compressed Sparse Row (CSR) or Compressed sparse columns (CSC) format. 
A variant of the variable Block Compressed Sparse Rows (or Columns) is used to store block-sparse matrices. 

SPARTA requires CUDA >=10.0 

# STRUCTURE

The files have the following structure

SPARTA
* include
* obj
* programs 
* src
* test   

each folder contains 
* general: files needed by all versions
* cuda: files needed by the cuda version
* mkl: files needed by the mkl version


# RUNNING A TEST

use `make` to create a test executable of the cuda test. The executable will be placed in programs/cuda. You can run it with different command line arguments to test different features.  
use 'source ./scripts/synthetic.sh' from the main folder to run and save some experiments. 

Options for the cuda_test:

* -i: select input example
* * 1: Random CSR
* * 3: Matrix Market (MTX) file
* * 4: Random Variable Block matrix
      
* -a: algorithm selection
* * -1: all
* * 1: cublas gemm
* * 2: VBSmm
* * 3: VBSmm with no zeros
* * 4: VBSmm with asymmetric hash-angle reordering
* * 5: cusparse spmm

* -b: density of blocks (% of nonzero blocks) (only for i = 4)

* -f: source file (only for i = 2, 3)

* -m: first matrix rows

* -n: second matrix columns

* -k: first matrix columns

* -p: size of VBS blocks

* -q: density of entries, in-block density (% of nonzero entries. if i = 4, % of nonzeros inside each nonzero block)

* -r: number of experiment repetitions

* -s: scramble input matrix. 
* * 0: no scramble. 
* * 1: scramble rows

* -S: random seed;

* -v: verbose level; ( -1: repeatead experiment format) 

* -w: warmup repetitions
