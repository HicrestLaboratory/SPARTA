# SPARTA
**SP**arse **A**ccele**R**ation on **T**ensor **A**rchitecture

The project aims to investigate new data structures and compression algorithms for exploting new architecture capabilities specifically designed for Deep Learning applications for **sparse and irregular computation**.

The repository contains code for compressing sparse matrices into dense block data-structure. 

Input sparse matrices are stored in Compressed Sparse Row (CSR) or Compressed sparse columns (CSC) format. 
Variable Block Compressed Sparse (Rows or Columns) is used to store block-sparse matrices. 

Running the code requires either Intel MKL library or the CUDA Toolkit.
MKL Library: https://software.intel.com/en-us/mkl

CUDA install;
Donwload the cuda toolkit (currently supported 10.0): https://developer.nvidia.com/cuda-downloads
Set the CUDA_PATH variable to your cuda path (default is /usr/local/cuda-10.0)


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

cuda_test:

* -i: select input example
* * 1: Random CSR
* * 2: SNAP Edgelist
* * 4: Random Variable Block matrix
      
* -a: algorithm selection
* * -1: all
* * 1: cublas gemm
* * 2: VBSmm
* * 3: VBSmm with no zeros
* * 4: VBSmm with asymmetric hash-angle reordering
* * 5: cusparse spmm

* -b: block sparsity (only for i = 4)

* -s: source file (only for i = 2, 3)

* -m: first matrix rows

* -n: second matrix columns

* -k: first matrix columns

* -p: size of VBS blocks

* -q: entries sparsity (if i = 4, block entries sparsity)

* -r: number of experiment repetitions

* -S: random seed;

* -v: verbose level; ( -1: repeatead experiment format) 

* -w: warmup repetitions
      
## RoadMAP
* Supporting Nvidia architecture (DONE!)
* New algorithm for Matrix reordering and compression to support no-coherent blocks (DONE!)
* New algorithm for sparse-sparse multiplication
**Stay Tuned!**
