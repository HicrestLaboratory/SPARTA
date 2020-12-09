# SPARTA
**SPAR**se acceleration on **T**ensor **A**rchitecture

The project aims to investigate new data structures and compression algorithms for exploting new architecture capabilities, specifically designed for deep learning, to accelerate **sparse and irregular** application, such as graph analytics.

The repository contains code for reordering and compressing sparse matrices into dense block data-structures.
The reordering algorithm matches rows (or columns) with similar patterns and attempts to build dense blocks. 
The similarity of patterns is first determined with a hash function, and then refined with a "angle algorithm", which matches patterns with high scalar product.

The repository also contains code for sparse-dense matrix-matrix multiplication that exploits the dense block data-structure.

Input sparse matrices are stored in Compressed Sparse Row (CSR) or Compressed sparse columns (CSC) format. 
A variant of the variable Block Compressed Sparse Rows (or Columns) is used to store block-sparse matrices. 

Running the code requires CUDA 10.0.

CUDA install
* Donwload the cuda toolkit (currently supported 10.0) and follow the instructions: https://developer.nvidia.com/cuda-downloads
* Set the CUDA_PATH variable to your cuda path (default is /usr/local/cuda-10.0)

# PRELIMINARY RESULTS
We have compared our routine with cusparse_spmm and cublas_gemm, the two main CUDA routines for sparse and dense matrix multiplication.

**CUSPARSE COMPARISON**
Preliminary results show that our routine is faster than cusparse_spmm when the density inside blocks is greater then around 2% (in this case, this corresponds to a total density of 0.2%)
![](/images/performance_experiment/VBS_vs_spmm_A8192_B_8192_fixed_blockdensity_0.1.jpg)


**CUBLAS COMPARISON**
Preliminary results show that our routine is faster than cublas_gemm when less than the 20% of blocks are nonzero. 
(both cublas and our routine treat nonzero blocks as dense, so changing the density inside blocks does not affect this result)
![](/images/performance_experiment_v2/VBS_vs_gemm_A4096_B_16384Block_size_128_varying_Block_density.jpg)


**PERFORMANCE LANDSCAPE**
The image below shows the faster algorithm for each data point when both the in-block density and the density of blocks vary.
For matrices in the green zone, SPARTA is the fastest choice. 
![](/images/performance_experiment_v2/scatter_performance_plot.jpg)

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

use `make test_cublas_VBS` to create a test executable of the cuda test. The executable will be placed in programs/cuda. You can run it with different command line arguments to test different features.  
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
* New algorithm for Matrix reordering and compression to support non-coherent blocks (DONE!)
* New algorithm for sparse-sparse multiplication
* Bring back support for MKL

**Stay Tuned!**
