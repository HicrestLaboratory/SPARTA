# SPARTA
**SP**arse **A**ccele**R**ation on **T**ensor **A**rchitecture

The repository contains some tools to manipulate CSR matrices, stored in SparMat objects. 
Sparse matrices can be read from Snap edgelist or MatrixMarket files into SparMat objects.

The make_sparse_blocks function rearranges the rows of a sparse matrix so that it can be separated in dense(r) blocks.
It also generates a variableSparseBlock (VBSparMat) matrix that stores these blocks efficiently.

Test.cpp shows an example of how reading, reorder and conversion work for a .mtx matrix.
It also generate some analysis of the new block structure. 
