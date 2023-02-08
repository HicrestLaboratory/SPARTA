#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <random> //std::shuffle
#include <stdexcept>

#include <algorithm>    // std::sort
#include <numeric> //std::iota
#include "matrices.h"
#include "utilities.h"

using namespace std;

void CSR::clean()
{
    /*----------------------------------------------------------------------
    | Free up memory allocated for CSR structs.
    |--------------------------------------------------------------------*/

    if (rows + cols <= 1) return;

    if (ma && !pattern_only) 
    {
        for(intT i = 0; i < rows; i++)
        {
            if (ma[i]) delete[] ma[i];
        }
    }

    if (ja) 
    {
        for(intT i = 0; i < rows; i++)
        {
            if (ja[i]) delete[] ja[i];
        }
    }
    
    if (nzcount) delete[] nzcount;
    
    rows = 0;
    cols = 0; 
}


void CSR::multiply(DataT* mat_B, intT B_cols, DataT_C* mat_C)
{
    //A*B, A is CSR, B and C column-wise storage

    for(intT i = 0; i < rows; i++)
    {
            for (intT nz = 0; nz < nzcount[i]; nz++)
            {
                intT col = ja[i][nz];
                DataT val = pattern_only?1:ma[i][nz];
                for(intT j = 0; j < B_cols; j++)
                {
                    mat_C[i + j*rows] += val*mat_B[col + j*rows];
                }
            }
    }
}

void CSR::permute_rows(vector<intT> permutation)
//permute the rows of the matrix according to "permutation"
{   
    if (permutation.size() != rows)
        throw std::invalid_argument("CSR.permute_rows argument bust have same lenght as rows");

    permute(ja,permutation);
    if (!pattern_only) permute(ma, permutation);
    permute(nzcount, permutation);
}

void CSR::reorder(vector<intT> grouping)
//permute the rows so that row i and row j are adjacent if grouping[i] == grouping[j]
{
    if (grouping.size() != rows)
        throw std::invalid_argument("CSR.reorder argument bust have same lenght as rows");

    vector<intT> v = get_permutation(grouping);
    permute_rows(v);
}

void CSR::reorder_by_degree(bool descending)
{
    auto asc_comparator = [&](int i, int j)
    {
        return nzcount[i] < nzcount[j];
    };

    auto desc_comparator = [&](int i, int j)
    {
        return nzcount[i] >= nzcount[j];
    };

    vector<intT> v(rows);
    iota(v.begin(), v.end(), 0);

    if (descending) sort (v.begin(), v.end(), desc_comparator);
    else sort (v.begin(), v.end(), asc_comparator);
    permute_rows(v);
}

void CSR::scramble()
{
    //randomly permute rows (TODO better randomness)
    vector<intT> v(rows);
    iota(v.begin(), v.end(), 0);

    random_shuffle(v.begin(),v.end());

    permute_rows(v);
}


void CSR::get_VBR_nzcount(const vector<intT> &grouping, intT block_col_size, intT& VBR_nzcount, intT& VBR_nzblocks_count, float& VBR_average_height)
{
    vector<intT> row_partition = get_partition(grouping);
    vector<intT> row_permutation = get_permutation(grouping);
    return get_VBR_nzcount(row_partition, row_permutation, block_col_size, VBR_nzcount, VBR_nzblocks_count, VBR_average_height);
}

void CSR::get_VBR_nzcount(const vector<intT> &row_partition, const vector<intT> &row_permutation, intT block_col_size, intT& VBR_nzcount, intT& VBR_nzblocks_count, float& VBR_average_height)
{
    
    intT block_cols = cols/block_col_size + 1;
    intT block_rows = row_partition.size() - 1;

    intT VBR_nzblocks_count = 0;
    intT total_blocks_height = 0;

    //copy data block_row by block_row
    for(intT ib = 0; ib < block_rows; ib++)
    {
        vector<bool> nonzero_flags(block_cols, false);
        intT row_block_size = row_partition[ib+1] - row_partition[ib];

        //flag nonzero blocks
        for (intT i_reordered = row_partition[ib]; i_reordered < row_partition[ib+1]; i_reordered++)
        {
            intT i = row_permutation[i_reordered];
            for (intT nz = 0; nz < nzcount[i]; nz++)
            {
                intT j = ja[i][nz];
                nonzero_flags[j/block_col_size] = true;
            }
        }

        for (intT jb = 0; jb < nonzero_flags.size(); jb++)
        {
            if (nonzero_flags[jb] == 1)
            {
                intT tmp_block_col_size = block_col_size - cols%block_col_size; //accounts for the last (possibly shorter) block
                VBR_nzcount += tmp_block_col_size*row_block_size;
                VBR_nzblocks_count++;
                total_blocks_height += row_block_size;
            }
        }
    }

    VBR_average_height = ((float) total_blocks_height)/VBR_nzblocks_count;
}

void CSR::read_from_edgelist(ifstream& infile, string delimiter, bool pattern_only)
//reads edgelist into the CSR.
{
    this->pattern_only = pattern_only;

    intT last_node = -1;
    intT current_node;
    string temp;
    vector<vector<intT>> pos_holder;
    vector<vector<DataT>> val_holder;
    intT max_column = 0;
    intT i = -1; 
    DataT val;
    intT total_nonzeros = 0;

    while (infile.peek() == '#' or infile.peek() == '%') infile.ignore(2048, '\n');
    while (getline(infile, temp)) {

        total_nonzeros++;
        int del_pos = temp.find(delimiter);
        int del_size = delimiter.length();

        string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);
        temp.erase(0, del_pos + del_size);
        
        del_pos = temp.find(delimiter);
        string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        intT child = stoi(second_node_string);
        max_column = std::max(max_column, child);
	    
        if (not pattern_only)
        {
            temp.erase(0, del_pos + del_size);
            del_pos = temp.find(delimiter);
            string val_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
            val = stof(val_string);
        }


        //fill with empty lines to reach current row
        if (current_node > i)
        {
            while (i < current_node)
            {
                vector<intT> new_pos_row;
		        pos_holder.push_back(new_pos_row);
                if (not pattern_only)
                {
                    vector<DataT> new_val_row;
                    val_holder.push_back(new_val_row);
                }
                i++;
            }
        }
        else if (current_node < i)
            throw std::invalid_argument("CSR.read_from_edgelist indices must be in ascending order");

        pos_holder[i].push_back(child);
    	if (not pattern_only) val_holder[i].push_back(val);
    }


    rows = pos_holder.size();
    cols = max_column + 1;
    nzcount = new intT[rows];
    ja = new intT*[rows];
    if (not pattern_only) 
    {
        ma = new DataT*[rows];
    }

    for (intT i = 0; i < pos_holder.size(); i++)
    {
        auto row_pos = pos_holder[i];
        nzcount[i] = row_pos.size();
        ja[i] = new intT[nzcount[i]]; 
        std::copy(row_pos.begin(), row_pos.end(), ja[i]);
	    pos_holder[i].clear();
	    
        if (not pattern_only)
        {
            auto row_val = val_holder[i];
            ma[i] = new DataT[nzcount[i]];
            std::copy(row_val.begin(), row_val.end(), ma[i]);
            val_holder[i].clear();
        }
    }

    pos_holder.clear();
    val_holder.clear();
}

void CSR::print(intT verbose)
{
    if (verbose > 0)
    {
        cout << "PRINTING A CSR MATRIX (arrays only)" << endl;
        cout << "ROWS: " << rows << " COLS: " << cols << " PATTERN_only: " << pattern_only << endl; 
        cout << "NZ: " << nztot() << endl;
    }
    if (verbose > 1)
    {
        cout << "JA:" << endl;
        for (intT i = 0; i < rows; i++)
        {
            cout << "-- ";
            for (intT j = 0; j < nzcount[i]; j++)
            {
                cout << ja[i][j] << " ";
            }
            cout << endl;

        }

        if (!pattern_only)
        {        
            cout << "MA:" << endl;
            for (intT i = 0; i < rows; i++)
            {
                cout << "-- ";
                for (intT j = 0; j < nzcount[i]; j++)
                {
                    cout << ma[i][j] << " ";
                }
                cout << endl;
            }
        }

        cout << "NZCOUNT:" << endl;
        for (intT i = 0; i < rows; i++)
        {
                cout << nzcount[i] << " ";
        }
        cout << endl;

        cout << "PRINTING A CSR MATRIX" << endl;
        cout << "ROWS:" << rows << " COLS:" << cols << " PATTERN_only:" << pattern_only << endl; 
        //loop through rows
        for (intT i = 0; i < rows; i++)
        {
            intT j = 0;
            for (intT nzs = 0; nzs < nzcount[i]; nzs++) 
            {
                intT nz_column = ja[i][nzs]; //find column (row) index of next nonzero element
                
                DataT elem;
                if (!pattern_only) elem = ma[i][nzs]; //value of that element;
                else elem = 1;
            
                while (j < nz_column)
                {
                    j++;
                    cout << 0 << " ";
                }
                cout << elem << " ";
                j++;
            }
            while (j < cols)
            {
                j++;
                cout << 0 << " ";
            }

        cout << endl;
        }
    }
}
