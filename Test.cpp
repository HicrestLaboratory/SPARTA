#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include "globheads.h"
#include "protos.h"
#include "utilities.h"


int main() {

	Mat mat;
	SparMat spmat;

	//create a random CSR matrix
	const int row_dimension = 1000;
	float sparsity = 0.0001;
	double eps = 0.1;

	random_sparse_mat(mat, row_dimension, sparsity); //generate random Mat (has issues with big numbers)
	convert_to_CSR(mat, spmat);

	//read from edgelist
	//read_snap_format(spmat, "testgraph.txt");
	//matprint(spmat);

	//convert_from_CSR(spmat, mat);

	int nBlock;
	int *nB = NULL, *perm = NULL;
	double *t_hash = NULL, *t_angle = NULL;
	VBSparMat vbmat;

	//find block decomposition
	//TODO only structure, no data (for unweighted graphs)
	if (init_blocks(&spmat, &nBlock, &nB, &perm, eps, t_hash, t_angle) != 0) {
		cout << "ERROR: COULD NOT CREATE PERMUTATION. Try with another epsilon" << endl;
		return -1;
	}
	//permute the sparse mat
	permute(spmat, perm);

	//transform spmat into variable block mat
	int ierr = csrvbsrC(1, nBlock, nB, &spmat, &vbmat);
	if (ierr != 0) cout << "error occurred while creating block matrix" << endl;

	ofstream CSV_out;
	CSV_out.open("output.txt");

	//TODO multiple entries
	string CSV_header = "MatrixSize,OriginalSparsity,Divisions,NonzeroBlocks,AvgBlockHeight,AvgBHError,AvgBlockLength,AvgBLError,NonzeroAreaFraction,AverageBlockPopulation,ABPError,NewSparsity";
	CSV_out << CSV_header << endl;
	features_to_CSV(&vbmat, CSV_out, true);
	CSV_out.close();

	//vbmatvec(&vbmat,x,y);
}
