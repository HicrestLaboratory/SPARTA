# SPARTA
**SPAR**se acceleration on **T**ensor **A**rchitecture

The project aims to investigate new data structures and compression algorithms for exploting new architecture capabilities, specifically designed for deep learning, to accelerate **sparse and irregular** application, such as graph analytics.

The repository contains code for reordering and compressing sparse matrices into dense block data-structures.
The reordering algorithm matches rows (or columns) with similar patterns and attempts to build dense blocks. 
The similarity of patterns is first determined with a hash function, and then refined with a "angle algorithm", which matches patterns with high scalar product.

The repository also contains code for sparse-dense matrix-matrix multiplication that exploits the dense block data-structure.

Input sparse matrices are stored in Compressed Sparse Row (CSR) or Compressed sparse columns (CSC) format. 
A variant of the variable Block Compressed Sparse Rows (or Columns) is used to store block-sparse matrices. 

Running the code requires either Intel MKL library (currently not supported) or the CUDA Toolkit.

MKL Library: 
https://software.intel.com/en-us/mkl

CUDA install
* Donwload the cuda toolkit (currently supported 10.0) and follow the instructions: https://developer.nvidia.com/cuda-downloads
* Set the CUDA_PATH variable to your cuda path (default is /usr/local/cuda-10.0)

The files have the following structure

SPARTA
* include
* obj
* programs 
* src   

each folder contains 
* general: files needed by all versions
* cuda: files needed by the cuda version
* mkl: files needed by the mkl version

running `make cuda_test` or `make mkl_test` (currently not supported) will create a test executable of the cuda version or the mkl version respectively. The executable will be placed in programs/cuda or programs/mkl. You can run it with different command line arguments to test different features.  

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

* -q: density of entries (% of nonzero entries. if i = 4, % of nonzeros inside each nonzero block)

* -r: number of experiment repetitions

* -s: scramble input matrix. 
* * 0: no scramble. 
* * 1: scramble rows

* -S: random seed;

* -v: verbose level; ( -1: repeatead experiment format) 

* -w: warmup repetitions
      
## RoadMAP
* Supporting Nvidia architecture (DONE!)
* New algorithm for Matrix reordering and compression to support no-coherent blocks (DONE!)
* New algorithm for sparse-sparse multiplication
* Bring back support for MKL

**Stay Tuned!**
