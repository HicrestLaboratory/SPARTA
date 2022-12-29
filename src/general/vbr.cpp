#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>    // std::sort, std::count
#include <numeric> //std::iota
#include "matrices.h"
#include "utilities.h"

using namespace std;

void VBR::clean()
{
    if (rows + cols <= 1) return;
    
    if (mab) delete[] mab;
    
    if (jab) delete[] jab;

    if (nzcount) delete[] nzcount;

    if (row_part) delete[] row_part;

    rows = 0;
    cols = 0;
    block_rows = 0;
    block_cols = 0;
    nztot = 0;
    block_col_size = 0;
}

DataT* VBR::get_block_start(intT row_block_idx)
{
    //returns the position in the mab vector where entries from a row block are stored. 
    DataT* ptr = mab;
    if (row_block_idx >= block_rows) 
    {
        cerr << " [RANGE ERROR] VBR::get_block_start" << endl;
        return mab;
    }

    for (intT ib = 0; ib < row_block_idx; ib++)
    {
        ptr += (row_part[ib+1] - row_part[ib])*nzcount[ib]*block_col_size;
    }
    return ptr;

}

void VBR::print(int verbose)
{
    intT* jab_ptr = jab;
    for (intT ib = 0; ib < block_rows; ib++)
    {   
        intT row_block_size = row_part[ib+1] - row_part[ib];
        DataT* data_pointer = get_block_start(ib);

        for (intT i = 0; i < row_block_size; i++)
        {
            
            intT jb = 0;
            for (intT nzb = 0; nzb < nzcount[ib]; nzb++)
            {

                intT nz_jb = jab_ptr[nzb]; //position of nonzero block
                while (jb < nz_jb)
                {
                    for (intT j = 0; j < block_col_size; j++) cout << "0 ";
                    cout << "| ";
                    jb++;
                }
                jb++;
                //print mab  slice
                for (intT j = 0; j < block_col_size; j++) 
                {
                    DataT d = data_pointer[nzb*block_col_size*row_block_size + i*block_col_size + j];
                    cout << d << " " ;
                }    
                cout << "| ";
            }

            while (jb < block_cols)
            {
                for (intT j = 0; j < block_col_size; j++) cout << "0 ";
                {
                    cout << "| ";
                }
                jb++;
            }
            cout << endl;
        }
        cout << endl;
        jab_ptr += nzcount[ib]; //move jab ptr to start of next block_row
    }
    cout << endl;
}


int VBR::partition_check(const vector<intT> &candidate_part)
{
//checks that a vector is a valid partition; anything different than 0 is an error;
    if (candidate_part.size() == 0) return 1;
    if (candidate_part.back() != rows) return 2;
    for (intT i = 1; i < candidate_part.size(); i++) 
    {
        if (candidate_part[i] < candidate_part[i - 1]) return 3;
    }
    return 0;
}

void VBR::fill_from_CSR_inplace(const CSR& cmat,const vector<intT> &grouping, intT block_size)
{
    //fill the VBR with entries from a CSR, with rows permuted and grouped according to grouping.

    vector<intT> row_partition = get_partition(grouping);
    vector<intT> row_permutation = get_permutation(grouping);

    rows = cmat.rows;
    cols = cmat.cols;
    block_col_size = block_size;
    
    block_cols = cols/block_size;
    block_rows = row_partition.size() - 1;

    //partition check

    if (int a = partition_check(row_partition) != 0)
    {
        cerr << "PARTITION CHECK ERROR: partition check failed with error " << a << endl;

    }
    row_part = new intT[row_partition.size()];
    copy(row_partition.begin(),row_partition.end(), row_part);

    nzcount = new intT[block_rows]{0};

    vector<DataT> mab_vec;
    vector<intT> jab_vec;

    //copy data block_row by block_row
    for(intT ib = 0; ib < block_rows; ib++)
    {
        vector<bool> nonzero_flags(block_cols, false);
        intT row_block_size = row_part[ib+1] - row_part[ib];

        //flag nonzero blocks
        for (intT i_reordered = row_part[ib]; i_reordered < row_part[ib+1]; i_reordered++)
        {
            intT i = row_permutation[i_reordered];
            for (intT nz = 0; nz < cmat.nzcount[i]; nz++)
            {
                intT j = cmat.ja[i][nz];
                nonzero_flags[j/block_col_size] = true;
            }
        }

        //fill jab vector
        for (intT jb = 0; jb < block_cols; jb++)
        {
            if (nonzero_flags[jb]) jab_vec.push_back(jb);
        }

        //fill nzcount
        nzcount[ib] = std::count(nonzero_flags.begin(), nonzero_flags.end(), true);


        //fill mab vector
        intT current_mab_size = mab_vec.size(); 
        mab_vec.resize(current_mab_size + nzcount[ib]*row_block_size*block_col_size, 0);
        for (intT i_reordered = row_part[ib]; i_reordered < row_part[ib+1]; i_reordered++)
        {
            intT i = row_permutation[i_reordered];
            for (intT nz = 0; nz < cmat.nzcount[i]; nz++)
            {
                intT j = cmat.ja[i][nz];
                DataT d;
                if (cmat.job == 0) d = 1;
                else d = cmat.ma[i][nz];

                //find position of d in mat
                intT j_block_position = j/block_size;
                intT tmp_block_count = std::count(nonzero_flags.begin(), nonzero_flags.begin() + j_block_position, true); //how many nz_blocks before current
                intT tmp_mab_pos = current_mab_size + tmp_block_count*block_col_size*row_block_size + (i_reordered - row_part[ib])*block_col_size + j%block_size; 
             
                mab_vec[tmp_mab_pos] = d;
            }
        }

    }

    nztot = mab_vec.size();
    jab = new intT[jab_vec.size()];
    mab = new DataT[mab_vec.size()];
    copy(jab_vec.begin(), jab_vec.end(),jab);
    copy(mab_vec.begin(),mab_vec.end(), mab);
}

void VBR::fill_from_CSR(const CSR& cmat,const vector<intT> &row_partition, intT block_size)
{
    //fill the VBR with entries from a CSR, with rows appearing in the same order and divided according to row_partition.

    rows = cmat.rows;
    cols = cmat.cols;
    block_col_size = block_size;
    
    block_cols = cols/block_size;
    block_rows = row_partition.size() - 1;


    //partition check

    if (int a = partition_check(row_partition) != 0)
    {
        cerr << "PARTITION CHECK ERROR: partition check failed with error " << a << endl;

    }
    row_part = new intT[row_partition.size()];
    copy(row_partition.begin(),row_partition.end(), row_part);

    nzcount = new intT[block_rows]{0};

    vector<DataT> mab_vec;
    vector<intT> jab_vec;

    //copy data block_row by block_row
    for(intT ib = 0; ib < block_rows; ib++)
    {
        vector<bool> nonzero_flags(block_cols, false);
        intT row_block_size = row_part[ib+1] - row_part[ib];

        //flag nonzero blocks
        for (intT i = row_part[ib]; i < row_part[ib+1]; i++)
        {
            for (intT nz = 0; nz < cmat.nzcount[i]; nz++)
            {
                intT j = cmat.ja[i][nz];
                nonzero_flags[j/block_col_size] = true;
            }
        }

        //fill jab vector
        for (intT jb = 0; jb < block_cols; jb++)
        {
            if (nonzero_flags[jb]) jab_vec.push_back(jb);
        }

        //fill nzcount
        nzcount[ib] = std::count(nonzero_flags.begin(), nonzero_flags.end(), true);


        //fill mab vector
        intT current_mab_size = mab_vec.size(); 
        mab_vec.resize(current_mab_size + nzcount[ib]*row_block_size*block_col_size, 0);
        for (intT i = row_part[ib]; i < row_part[ib+1]; i++)
        {
            for (intT nz = 0; nz < cmat.nzcount[i]; nz++)
            {
                intT j = cmat.ja[i][nz];
                DataT d;
                if (cmat.job == 0) d = 1;
                else d = cmat.ma[i][nz];

                //find position of d in mat
                intT j_block_position = j/block_size;
                intT tmp_block_count = std::count(nonzero_flags.begin(), nonzero_flags.begin() + j_block_position, true); //how many nz_blocks before current
                intT tmp_mab_pos = current_mab_size + tmp_block_count*block_col_size*row_block_size + (i - row_part[ib])*block_col_size + j%block_size; 
             
                mab_vec[tmp_mab_pos] = d;
            }
        }

    }

    nztot = mab_vec.size();
    jab = new intT[jab_vec.size()];
    mab = new DataT[mab_vec.size()];
    copy(jab_vec.begin(), jab_vec.end(),jab);
    copy(mab_vec.begin(),mab_vec.end(), mab);
    cout << "mab: ";
    print_vec(mab_vec);
    cout << "jab: ";
    print_vec(jab_vec);
}
