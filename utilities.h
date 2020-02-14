#pragma once

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <random>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <map>
#include <set>
#include <iomanip> 


#include "globheads.h"
#include "protos.h"
#include "mkl.h"

using namespace std;

typedef map<int, set<int> > Graphmap;


struct Mat {
	int row_size;
	vector<double> vec;
};

void random_sparse_mat(Mat & mat, int N, double sparsity);

void graph_print(Graphmap & gmap);

void read_snap_format(Graphmap & gmap, string filename);

int isProper(const Graphmap & gmap, bool mirroring);

void MakeUndirected(Graphmap & gmap);

void MakeProper(Graphmap & gmap);

void write_snap_format(Graphmap & gmap, string filename);

void convert_from_CSR(SparMat & spmt, Mat & mat);

void fill_CSR(SparMat & spmt, int n, const vector<int>& vec_nzcount, const vector<int>& vec_ja, const vector<double>& vec_ma);

void convert_to_CSR(const Mat & mat, SparMat & spmt);

void convert_to_CSR(const Graphmap & gmap, SparMat & spmt);

void permute(SparMat & spmt, int * perm);

void matprint(const SparMat & spmt);

void matprint(const Mat & mat);

void matprint(const double* mat, const int n,const int m);

void matprint(const VBSparMat &vbmat);

void convert_to_MKL(SparMat &spmt, sparse_matrix_t &A);

void read_snap_format(SparMat & spmt, string infilename);

void read_mtx_format(SparMat & spmt, string infilename);

void extract_features(const vbsptr vbmat, int & Msize, int & Bnum, vector<int>& BLsize, vector<int>& BHsize, vector<int>& Bsparsity);

void features_to_CSV(vbsptr vbmat, ofstream & CSV_out, int verbose);

int make_sparse_blocks(SparMat &spmt, VBSparMat &vbmat,double eps);

void convert_to_col_major(double *X, double *Y, const int n, const int m);

void convert_to_row_major(double *X, double *Y, const int n, const int m);

bool are_equal(const double *X,const double* Y,const int m, const double eps = 0.0);

void block_mat_multiply(const VBSparMat &VBMat, double *X, const int k, double *Y);

//-----------------------------------------------------------------------ARRAY AND VECTOR UTILITIES----------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------

void randomvec(vector<double>& v, int n);


//mean of a vector
template <class myType>
double mean(const vector<myType>& v) {
	return accumulate(v.begin(), v.end(), 0.) / v.size();
}
//-----------------------------------------------------------------//


//Standard deviation of a vector
template <class myType>
double stdv(const vector<myType>& v) {

	if (v.size() < 2) { return 0; }

	double u = mean(v);
    auto lambda = [u](double a, myType b) {
		return a + pow(((double)b - u), 2);
	};

	double a = accumulate(v.begin(), v.end(), 0., lambda);
	return sqrt(a / (v.size() - 1));
}
//-----------------------------------------------------------------//


//permutes an array of n pointers (original) according to a permutation (perm);
template <class myType>
void permute(myType* original, int* perm, int n) {
	int i = 0;
	myType temp;
	int cur_row, next_row;

	//Do the permutations cycle by cycle, until all elements have been permuted.
	while (i < n && perm[i] >= 0) {//LI: all elements up to i appeared in a cycle
		temp = original[i];
		cur_row = i;
		next_row = perm[cur_row];

		//one cycle of permutations
		while (perm[cur_row] >= 0) {//check if we are back to the start of the cycle
			next_row = perm[cur_row];
			perm[cur_row] = (perm[cur_row] + 1)*(-1);//flag perm element as executed
			original[cur_row] = original[next_row];
			cur_row = next_row;
		}
		//end cycle and procede to next element
		original[cur_row] = temp;
		i++;
	}

	//reconstruct perm (remove flags)
	for (i = 0; i < n; i++) {
		perm[i] = -perm[i] - 1;
	}


}

//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------


