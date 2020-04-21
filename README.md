# SPARTA
**SP**arse **A**ccele**R**ation on **T**ensor **A**rchitecture

The project aims at investigating new data structures and compression algorithms for exploting new architecture capabilities specifically designed for Deep Learning applications for **sparse and irregular computation**.

The repository contains code for compressing sparse matrices into dense block data-structure. 


Input sparse matrices are stored in Compressed Sparse Row (CSR) format. 
Variable Block Compressed Sparse Row is used to store block-sparse matrices. 

Running the code requires either Intel MKL library or the CUDA Toolkit.
MKL Library: https://software.intel.com/en-us/mkl
CUDA: https://developer.nvidia.com/cuda-downloads

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

running `make cuda_test` or `make mkl_test` will create a test executable of the cuda version or the mkl version respectively. The executable will be placed in programs/cuda or programs/mkl. You can run it with different command line arguments to test different features.  

* `-i` : select type of input
* *  `1`: Random CSR
* *  `2`: SNAP Edgelist
* *  `3`: MTX Format
* *  `4`: Random Variable Block matrix
* `-s`: provide custom source file (has only effect for example 2 and 3)
* `-k`: input matrix sparsity (has only effect for example 1 and 4)
* `-n`: input matrix dimension (has only effect for example 1 and 4)
* `-o`: number of column of output matrix
* `-e`: epsilon used for marix reordering (must be in [0,1]. Smaller epsilon will make larger blocks);  


## RoadMAP
* Supporting Nvidia Tensor Core
* New algorithm for Matrix Compression to support no-coherent blocks
* Graph Analytics applications
**Stay Tuned!**
