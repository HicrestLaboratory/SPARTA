#pragma once

#include "matrices"

using namespace std;

void CSR::reorder(vector<intT> permutation)
{
}

void CSR::read_from_edgelist(ifstream& infile, string delimiter, bool pattern_only)
{
    intT last_node = -1;
    intT current_node;
    std::string temp;
    intT max_column = 0;
    vector<intT> pos_holder;
    vector<DataT> val_holder;
    intT i = -1; 
    
    while (infile.peek() == '#' or infile.peek() == '%') infile.ignore(2048, '\n');
    //while (new node equal last):
    while (getline(infile, temp)) {

        int del_pos = temp.find(delimiter);
        int del_size = delimiter.length();

        string first_node_string = temp.substr(0, del_pos); //retrieve the part of the string before the delimiter
        current_node = stoi(first_node_string);
        temp.erase(0, del_pos + del_size);
        
        del_pos = temp.find(delimiter);
        std::string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        intT child = stoi(second_node_string);
        max_column = std::max(max_column, child);
        temp.erase(0, del_pos + del_size);

	if (!pattern_only) 
	{
		del_pos = temp.find(delimiter);
        	std::string second_node_string = temp.substr(0, del_pos); //retrieve the part of the string after the delimiter
        	DataT weigth = stof(second_node_string);
	}


	//if current_node != last_node:
	//	save

        if (current_node > i)
        {
            while (i < current_node)
            {
                i++;
            }
        }
        else if (current_node < i)
        {
            std::cerr << "CANNOT READ MATRIX. INDICES MUST INCREASE" << std::endl;
            return 1;
        }
        holder[i].push_back(child);
    }

    cmat.fmt = cmat_fmt;
    intT main_dim = holder.size();
    intT second_dim = max_column + 1;
    cmat.rows = cmat_fmt ? second_dim : main_dim;
    cmat.cols = cmat_fmt ? main_dim : second_dim;
    cmat.nzcount = new intT[main_dim];
    cmat.ja = new intT * [main_dim];
    cmat.ma = new DataT * [main_dim];

    for (intT i = 0; i < holder.size(); i++)
    {
        auto row = holder[i];
        cmat.nzcount[i] = row.size();
        cmat.ja[i] = new intT[row.size()];
        cmat.ma[i] = new DataT[row.size()];

        std::copy(row.begin(), row.end(), cmat.ja[i]);
        std::vector<DataT> temp_vec(row.size(), 1.);
        std::copy(temp_vec.begin(), temp_vec.end(), cmat.ma[i]); //entries = 1s. unweigthed
    }

    holder.clear();
    return 0;
}

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