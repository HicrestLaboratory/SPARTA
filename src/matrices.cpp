#pragma once

#include "matrices"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

void CSR::print(ofstream& outfile)
{
    //loop through rows
    for (intT i = 0; intT i < rows; i++)
    {
	intT last_col = 0;
        for (intT nzs = 0; nzs < nzcount[i]; nzs++) 
        {
	     
            nz_column = ja[i][nzs]; //find column (row) index of next nonzero element
            DataT elem = 1;
	    if (job == 1) elem = ma[i][nzs]; //value of that element;
	    for (intT j = last_col; j < nz_column; j++)
	    {
		    outfile << 0 << " ";
	    }
	    outfile << elem << " ";
	    last_col = nz_column;
        }
	
    	for (intT j = last_col; j < cols; j++)
    	{
	    outfile << 0 << " ";
    	}	
	outfile << endl;
    }

}

void CSR::reorder(vector<intT> permutation)
{
}

void CSR::read_from_edgelist(ifstream& infile, string delimiter = "\t", bool pattern_only = true)
{

    intT last_node = -1;
    intT current_node;
    string temp;
    vector<vector<intT>> pos_holder;
    vector<vector<DataT>> val_holder;
    intT max_column = 0;
    intT i = -1; 
    
    while (infile.peek() == '#' or infile.peek() == '%') infile.ignore(2048, '\n');
    while (getline(infile, temp)) {

        int del_pos = temp.find(delimiter);
        int del_size = delimiter.length();

        string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);
        temp.erase(0, del_pos + del_size);
        
        del_pos = temp.find(delimiter);
        string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        intT child = stoi(second_node_string);
        max_column = std::max(max_column, child);
	    
	del_pos = temp.find(delimiter);
        string val_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        DataT val = stof(val_string);
	    
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
        {
            std::cerr << "CANNOT READ MATRIX. INDICES MUST INCREASE" << std::endl;
            return 1;
        }
        pos_holder[i].push_back(child);
	val_holder[i].push_back(val);

    }

    rows = holder.size();
    cols = max_column + 1;
    nzcount = new intT[main_dim];
    ja = new intT * [main_dim];
    if (not pattern_only) ma = new DataT * [main_dim];

    for (intT i = 0; i < pos_holder.size(); i++)
    {
        auto row_pos = pos_holder[i];
        nzcount[i] = row_pos.size();
        ja[i] = new intT[row_pos.size()];
        std::copy(row_pos.begin(), row_pos.end(), ja[i]);
	pos_holder[i].clear();
	    
	    
	if (not pattern_only)
	{
	auto row_val = val_holder[i];
	ma[i] = new DataT[row_val.size()];
        std::copy(row_val.begin(), row_val.end(), ma[i]);
	val_holder[i].clear();
	}
    }

    pos_holder.clear();
    val_holder.clear();
    return 0;
}

void CSR::clean()
{
    /*----------------------------------------------------------------------
    | Free up memory allocated for CSR structs.
    |--------------------------------------------------------------------*/

    if (rows + cols <= 1) return;

    for (intT i = 0; i < rows; i++) {
        if (nzcount[i] > 0) {
            if (job) 
		{
		if (ma) delete[] ma[i];
		}
            delete[] ja[i];
        }
    }
    if (ma) delete[] ma;
    delete[] ja;
    delete[] nzcount;
}
