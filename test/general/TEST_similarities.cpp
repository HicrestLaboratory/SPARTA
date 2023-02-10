#include <fstream>
#include <iostream>
#include <algorithm>
#include "utilities.h"
#include "blocking.h"

using namespace std;

int main(int argc, char* argv[])
{
    CLineReader cli(argc, argv);
    intT block_size = cli.col_block_size_;


    vector<vector<intT>> tests;
    tests.push_back(vector<intT>{1,3}); 
    tests.push_back(vector<intT>{0,5,10}); 
    tests.push_back(vector<intT>{1,2,5,10,12,20}); 
    tests.push_back(vector<intT>{0,3,7,11}); 



    
    for (auto it = tests.begin(); it != tests.end(); it++)
        for (auto it2 = tests.begin(); it2 != tests.end(); it2++)
        {
            cout << "*******************************************************************" << endl;
            intT size_B = it2->size();
            intT row_B[size_B];
            copy(it2->begin(), it2->end(), row_B);

            float hamming_sim = HammingDistanceGroup(*it,1,row_B,size_B,1,block_size);
            float jaccard_sim = JaccardDistanceGroup(*it,1,row_B,size_B,1,block_size);

            float hamming_sim_openmp = HammingDistanceGroupOPENMP(*it,1,row_B,size_B,1,block_size);
            float jaccard_sim_openmp = JaccardDistanceGroupOPENMP(*it,1,row_B,size_B,1,block_size);

            cout << "block_size: " << block_size << endl;

            cout << "row A:  ";
            print_vec(*it);

            cout << "row B:  ";
            print_mat(row_B, 1, size_B, 1);

            cout << " Hamming: " << hamming_sim << endl;
            cout << " Hamming(OPENMP): " << hamming_sim_openmp << endl;
            
            if (!(hamming_sim == hamming_sim_openmp)) cout << "HAMMING TEST FAILED!" << endl;


            cout << " Jaccard: " << jaccard_sim << endl;
            cout << " Jaccard(OPENMP): " << jaccard_sim_openmp << endl;
            
            if (!(jaccard_sim == jaccard_sim_openmp)) cout << "JACCARD TEST FAILED!" << endl;



        }
}