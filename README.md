# SPARTA
**SP**arse **A**ccele**R**ation on **T**ensor **A**rchitecture

The repository contains code for compressing sparse matrices into dense block data-structure. This allow to manipulate (e.g. multiplicate) these matrices by esploiting on AI-specific architecture (e.g., Nvidia tensor core). 

Sparse matrices are efficiently stored in SparMat objects, which follow the Compressed Sparse Row (CSR) format. block-sparse matrices are strored in VBSparMat object which implement Variable Block Compressed Sparse Row. 

Currently, the code requires Intel MKL library.  Remember to configure the link command through
[Intel MKL Link Line Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/)

running `make` will create a test executable (from Test.cpp) in the programs/ folder. You can run it with different command line arguments to test different features.  

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


Source files are stored in the src folder:
* Test.cpp shows an example of how reading, conversion, rows reorder and sparse-dense multiplication work for several input matrix format.
* setblks.cpp contains tools to handle interactions (conversion, reorder) between SparMat and VBSparMat objects.
* utilities.cpp contains other tools (mainly interaction with SparMat objects). 

## RoadMAP
* Supporting Nvidia Tensor Core
* New algorithm for Matrix Compression to support no-coherent blocks
* Graph Analytics applications
**Stay Tuned!**
