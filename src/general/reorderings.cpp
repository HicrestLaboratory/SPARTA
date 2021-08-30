#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream> //ifstream

#include <cmath> // std::abs
#include <string>
#include <map>
#include <set>

#include <random>
#include <vector>
#include <algorithm>    // std::random_shuffle

#include "sparse_utilities.h"
#include "comp_mats.h"
#include "reorderings.h"
#include "input.h"

intT count_groups(intT* grp, intT grp_len)
{
    intT* perm = new intT[grp_len];
    sort_permutation(perm, grp, grp_len);
    intT last_grp = -1;
    intT groups = 0;

    for (intT idx = 0; idx < grp_len; idx++)
    {
        intT i = perm[idx];
        if (grp[i] != last_grp)
        {
            groups += 1;
            last_grp = grp[i];
        }
    }

    delete[] perm;
    return groups;

}

int grp_to_partition(intT* grp, intT grp_len, intT* partition)
{
    // IN: 
    //    grp: an array with group assignments
    //    grp_len: lenght of grp
    // OUT: 
    //    partition: partition similar entries of grp together (sorted by entry)
    intT* perm = new intT[grp_len];
    sort_permutation(perm, grp, grp_len);
    intT last_grp = -1;
    intT group = 0;

    for (intT idx = 0; idx < grp_len; idx++)
    {
        intT i = perm[idx];
        if (grp[i] != last_grp)
        {
            partition[group] = idx;
            group += 1;
            last_grp = grp[i];
        }

    }
    partition[group] = grp_len;

    delete[] perm;
}

int hash_permute(CSR& cmat, intT* compressed_dim_partition, intT* perm, intT* group, int mode)
{
    //finds a group structure and a permutation for the main dimension of a CSR mat
    //NOTE: this does NOT automatically permute the CSR

    // IN:
    //  cmat:        a matrix in CSR form
    //  block_size:  the number of elements to consider together when determining sparsity structure
    //              e.g if block_size = 8, every 8 element of secondary dimension will be considered nonzero if any of that is so
    //  mode:        0: at most one element per block is considered
    //              1: all elements in a block contribute to the hash
    // OUT:
    //  perm:        an array of length equal to cmat main dimension; stores a permutation of such dimension that leaves groups together
    //  group:       an array of length equal to cmat main dimension; for each main row, stores the row group it belongs to

    intT main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    intT* hashes = new intT[main_dim]; //will store hash values. The hash of a row (col) is the sum of the indices (mod block_size) of its nonzero entries

    for (intT i = 0; i < main_dim; i++)
    {
        group[i] = -1;

        hashes[i] = hash(cmat.ja[i], cmat.nzcount[i], compressed_dim_partition, mode); //calculate hash value for each row
    }

    sort_permutation(perm, hashes, main_dim); //find a permutation that sorts hashes


    intT* ja_0, * ja_1;
    intT len_0, len_1;
    intT tmp_group = -1;

    for (intT idx = 0; idx < main_dim; idx++) //scan main dimension in perm order and assign rows with same pattern to same group;
    {
        intT i = perm[idx]; //counter i refers to original order. Counter idx to permuted one. 

        if (group[i] == -1) //if row is still unassigned
        {

            tmp_group++; //create new group
            group[i] = tmp_group; //assign row to group

            ja_0 = cmat.ja[i]; //the row in compressed sparse format
            len_0 = cmat.nzcount[i];//the row length
            for (intT jdx = idx + 1; jdx < main_dim; jdx++)
            {
                intT j = perm[jdx]; //counter j refers to original order. Counter jdx to permuted one. 
                if (hashes[j] != hashes[i]) break; //only row with the same hashes must be compared, and they are consecutive in the perm order
                if (group[j] != -1) continue; //only unassigned row must be compared

                ja_1 = cmat.ja[j];
                len_1 = cmat.nzcount[j];

                if (check_same_pattern(ja_0, len_0, ja_1, len_1, compressed_dim_partition, mode))
                {
                    group[j] = tmp_group; //assign row j to the tmp_group
                }
            }
        }
    }

    delete[] hashes;
    sort_permutation(perm, group, main_dim); //stores in perm a permutation that sorts rows by group
}


intT hash(intT* arr, intT a_len, intT* block_partition, int mode)
{
    //evaluate hash function for a arbitrarily partitioned array of indices
    /* IN:
            arr: the array of indices (a compressed row)
            a_len: length of the array;
            block_partition: start position of block i; elements in the same block give the same contribution to hash
            mode:  0: at most one element per block contribute to hash
                   1: all elements in a block contribute to hash
        OUT:
            intT hash : the hash is the sum of the indices of nonzero blocks (indices are counted from 1, to avoid ignoring the 0 idx);
            */


    intT nzs = 0;
    intT hash = 0;

    intT prev_idx = -1;
    intT block_idx = 0;

    while (nzs < a_len)
    {
        while (arr[nzs] >= block_partition[block_idx + 1])
        {
            block_idx++;
        };

        nzs++;
        if ((block_idx == prev_idx) and (mode == 0)) //if mode is 0, only one element per block is considered in the hash sum;
        {
            continue;
        }

        hash += block_idx + 1;
        prev_idx = block_idx;
    }
    return hash;
}


int check_same_pattern(intT* arr0, intT len0, intT* arr1, intT len1, intT block_size, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DEFINED FOR BLOCKS OF FIXED DIMENSION block_size


    intT i = 0;
    intT j = 0;
    intT b_idx0 = 0;
    intT b_idx1 = 0;

    while ((i < len0) and (j < len1))
    {
        b_idx0 = arr0[i] / block_size;
        b_idx1 = arr1[j] / block_size;
        if (b_idx0 != b_idx1)
        {
            return 0;
        }

        i++;
        j++;

        if (mode == 0) //if mode=0, skip all entries in a block after the first one;
        {
            while ((i < len0) and (b_idx0 == arr0[i] / block_size))
            {
                i++;
            }
            while ((j < len1) and (b_idx1 == arr1[j] / block_size))
            {
                j++;
            }
        }

    }
    if ((i < len0) or (j < len1))
    {
        return 0;
    }
    else
    {
        return 1;
    }

}


int check_same_pattern(intT* arr0, intT len0, intT* arr1, intT len1, intT* block_partition, int mode)
{
    //check if two arrays of indices have the same pattern
    //PATTERN IS DIVIDED IN BLOCKS OF VARIABLE DIMENSION, GIVEN BY block_partition; 

    intT i = 0;
    intT j = 0;

    intT block_idx_0 = 0; //block_idx for arr0
    intT block_idx_1 = 0; //block_idx for arr1


    while ((i < len0) and (j < len1))
    {
        while (arr0[i] >= block_partition[block_idx_0 + 1])
        {
            block_idx_0++;
        }
        while (arr1[j] >= block_partition[block_idx_1 + 1])
        {
            block_idx_1++;
        }
        if (block_idx_0 != block_idx_1)
        {
            return 0;
        }

        i++;
        j++;

        if (mode == 0) //if mode=0, skip all entries in a block after the first one;
        {
            while ((i < len0) and (arr0[i] < block_partition[block_idx_0 + 1]))
            {
                i++;
            }
            while ((j < len1) and (arr1[j] < block_partition[block_idx_1 + 1]))
            {
                j++;
            }
        }

    }
    if ((i < len0) or (j < len1))
    {
        return 0;
    }
    else
    {
        return 1;
    }

}


int get_pattern(intT* arr0, intT len0, intT* block_partition, intT* pattern, int mode)
{
    /*get the pattern of a compressed array w.r.t. a given block_partition
    IN:
        arr0: the compressed array, storing index of nonzero elements of the uncompressed array.
        len0: length of the compressed array, i.e. # of nonzero elements in the uncrompressed array;
        block_partition, bp: the partition. Block i has boundary [ bp[i], bp[i+1] ) ;
        mode:      0: at most one element per block contributes
                   1: all elements in a block contribute

    OUT:
        pattern: the pattern (has length of block partition - 1)
    */

    intT i = 0;
    intT block_idx = 0; //block_idx for pattern
    intT in_block = 0;

    while (i < len0)
    {
        while (arr0[i] >= block_partition[block_idx + 1]) //check if idx out of block; in that case, procede to next block.
        {
            pattern[block_idx] = in_block;
            in_block = 0;
            block_idx++;
        }

        in_block += 1;
        if ((mode == 0) and (in_block > 1))
        {
            in_block = 1;
        }

        pattern[block_idx] = in_block;
        i++;

    }
    return 0;
}


int angle_method(CSR& cmat, float eps, intT* compressed_dim_partition, intT nB, intT* in_perm, intT* in_group, intT* out_group, int mode)
{
    /*
    COMPUTES A GROUPING AND PERMUTATION FOR A CSR MATRIX (ALONG ITS MAIN DIMENSION)
    GIVEN A STARTING PERMUTATION AND A GROUPING; USES ANGLE METHOD WITH A FIXED PARTITION OF THE COMPRESSED DIMENSION.

    IN:
        cmat: the CSR matrix
        compressed_dim_partition: the partition of the matrix along its compressed dimension. Element i is the start position of block i; last element is the length of the compressed dimension;
        nB: the number of blocks in the partition
        in_perm: a permutation of the main dimension (rows in the same group should be adjacent in this permutation);
        in_group: a grouping along the main dimension;
        eps: rows with cosine greater than eps will be merged.
        mode: 0: consider at most one element per block when evaluating sparsity patterns
              1: consider all the elements in a block (but not their order) when evaluating sparsity patterns

    OUT: out_group will store a new grouping (compatible with in_group)
    */

    intT main_dim = cmat.fmt == 0 ? cmat.rows : cmat.cols;
    intT second_dim = cmat.fmt == 0 ? cmat.cols : cmat.rows;

    for (intT i = 0; i < main_dim; i++)//initialize out_group
    {
        out_group[i] = -1;
    }

    intT idx, jdx;
    intT this_group = -1; //the (out_)group the current row is in. 
    intT in_this_grp;
    intT* this_pattern = new intT[nB]{ 0 }; //the pattern of the current row or group of rows

    intT that_group;//the (in_)group the compared row is in
    intT in_that_grp;
    intT* that_pattern = new intT[nB]{ 0 };

    intT i, j;

    for (intT idx = 0; idx < main_dim; idx++) //Loop through (groups of) rows. Each one is confronted with all the unpaired ones to find those that will be merged.
    {
        i = in_perm[idx];           //idx counts in the permuted order. i counts in the original order;

        if (out_group[i] == -1)     //only consider still ungrouped rows;
        {
            this_group++;               //create new group

            intT* arr0 = cmat.ja[i];
            intT len0 = cmat.nzcount[i];

            jdx = idx; //idx to compare the row with. will only compare with elements after current row;

            in_this_grp = 0;
            while ((jdx < main_dim) and (in_group[in_perm[jdx]] == in_group[i])) //check for elements in the same in_group of i (they are consecutive in in_perm)
            {
                out_group[in_perm[jdx]] = this_group;   //assign elements in the same in_group to the same out_group;
                jdx++;
                in_this_grp++; //keep count of the elements in the group
            }

            get_pattern(arr0, len0, compressed_dim_partition, this_pattern, mode); //get the row pattern (stores into this_pattern)




            if (mode == 1)
            {
                for (intT b = 0; b < nB; b++)
                {
                    this_pattern[b] *= in_this_grp; //multiply entries of the pattern with number of rows with the same pattern (same in_group)
                }
            }
            intT norm_0 = norm2(this_pattern, nB); //squared norm of the pattern

            while (jdx < main_dim) //loop through not-analyzed rows to be paired with the group
            {
                j = in_perm[jdx];
                that_group = in_group[j];
                bool merge = false;

                if (out_group[j] == -1)         //only consider ungrouped rows;
                {
                    intT* arr1 = cmat.ja[j];
                    intT len1 = cmat.nzcount[j];
                    that_group = in_group[j];

                    get_pattern(arr1, len1, compressed_dim_partition, that_pattern, mode); //get the row pattern (store into that_pattern)

                    intT norm_1 = norm2(that_pattern, nB); //get norm of the pattern

                    float scal = scalar_product(this_pattern, nB, that_pattern);
                    if (scal * scal > eps * norm_0 * norm_1) //if cosine is > than epsilon, allow merge of groups
                    {
                        merge = true;
                    }
                }

                in_that_grp = 0;
                while ((jdx < main_dim) and (in_group[in_perm[jdx]] == that_group)) //iterate over elements in the same group of the analyzed row j (j included)
                {
                    in_that_grp++;
                    if (merge) out_group[in_perm[jdx]] = this_group; //if merge was decided, add group to current group 
                    jdx++;
                }

                if (merge)
                {
                    for (intT b = 0; b < nB; b++)
                    {
                        if (mode == 0)
                        {
                            this_pattern[b] = std::max(this_pattern[b], that_pattern[b]);
                        }
                        else
                        {
                            this_pattern[b] += in_that_grp * that_pattern[b]; //update pattern with entries of merged rows
                        }
                    }
                    norm_0 = norm2(this_pattern, nB); //update norm of this pattern
                }
            }
        }
    }

    delete[] this_pattern;
    delete[] that_pattern;


}


int angle_hash_method(CSR& cmat, float eps, intT* compressed_dim_partition, intT nB, VBS& vbmat, int vbmat_blocks_fmt, int vbmat_entries_fmt, int mode)
{
    //create a VBS reordering the main (uncompressed) dimension of a CSR matrix according to the angle+hash algorithm
    //do not change the original array

    intT rows = cmat.rows;
    intT cols = cmat.cols;
    intT main_dim = (cmat.fmt == 0) ? rows : cols;


    intT* hash_perm = new intT[main_dim];
    intT* hash_grp = new intT[main_dim];

    hash_permute(cmat, compressed_dim_partition, hash_perm, hash_grp, mode);

    intT* angle_perm = new intT[main_dim];
    intT* angle_grp = new intT[main_dim];

    angle_method(cmat, eps, compressed_dim_partition, nB, hash_perm, hash_grp, angle_grp, mode);

    sort_permutation(angle_perm, angle_grp, main_dim); //find a permutation that sorts groups

    intT angle_main_grps;
    angle_main_grps = count_groups(angle_grp, main_dim);

    intT* angle_main_part = new intT[angle_main_grps + 1];

    grp_to_partition(angle_grp, main_dim, angle_main_part);


    CSR cmat_cpy;
    copy(cmat, cmat_cpy);


    permute_CSR(cmat_cpy, angle_perm, cmat_cpy.fmt); //permute the tmp CSR

    intT* row_part;
    intT row_blocks;
    intT* col_part;
    intT col_blocks;

    if (cmat.fmt == 0)
    {
        row_part = angle_main_part;
        row_blocks = angle_main_grps;
        col_part = compressed_dim_partition;
        col_blocks = nB;
    }
    else
    {
        col_part = angle_main_part;
        col_blocks = angle_main_grps;
        row_part = compressed_dim_partition;
        row_blocks = nB;
    }

    convert_to_VBS(cmat_cpy,
        vbmat,
        angle_main_grps, angle_main_part,
        col_blocks, col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt);

    //cleaning
    cleanCSR(cmat_cpy);
    delete[] hash_perm;
    delete[] hash_grp;
    delete[] angle_perm;
    delete[] angle_grp;
}

intT row_hash(intT* cols, intT len)
{
    if (len == 0) return -1;
    intT hash = 0;
    for (int i = 0; i < len; i++)
    {
        hash += cols[i];
    }
    return hash;
}

intT row_block_hash(intT* cols, intT len, intT block_size)
{
    if (len == 0) return -1;
    intT hash = 0;

    intT i = 0;
    intT block_idx;
    while (i < len)
    {
        block_idx = cols[i] / block_size;
        hash += block_idx;
        while (i < len && cols[i] / block_size == block_idx) i++;
    }

    return hash;
}


bool equal_rows(intT* cols_A, intT len_A, intT* cols_B, intT len_B)
{
    if (len_A != len_B) return false;
    else
    {
        intT i = 0;
        while (cols_A[i] == cols_B[i] && i < len_A)
        {
            i++;
        }
        if (i == len_A) return true;
        else return false;
    }
}

int hash_reordering(CSR& cmat, intT* groups, input_parameters &params)
{
    intT* hashes = new intT[cmat.rows]{ 0 };

    if (params.reorder_algo == "saad")
    {
        for (int i = 0; i < cmat.rows; i++)
        {
            hashes[i] = row_hash(cmat.ja[i], cmat.nzcount[i]);
        }
    }
    else if (params.reorder_algo == "saad_blocks")
    {
        for (int i = 0; i < cmat.rows; i++)
        {
            hashes[i] = row_block_hash(cmat.ja[i], cmat.nzcount[i], params.algo_block_size);
        }
    }

    intT* perm = new intT[cmat.rows]{ 0 };
    sort_permutation(perm, hashes, cmat.rows);

    intT current_group = 0;
    groups[perm[0]] = current_group;

    for (int ip = 1; ip < cmat.rows; ip++)
    {
        intT curr = perm[ip];
        intT prev = perm[ip - 1];
        if (hashes[curr] != hashes[prev])
        {
            current_group++;
        }
        else
        {
            if (!equal_rows(cmat.ja[curr], cmat.nzcount[curr], cmat.ja[prev], cmat.nzcount[prev])) current_group++;
        }
        groups[curr] = current_group;
    }
}


int assign_group(intT* in_group, intT* out_group, intT* perm, intT len, intT jp, intT new_group_idx)
{
    intT current_in_group = in_group[perm[jp]];
    while (jp < len && in_group[perm[jp]] == current_in_group)
        {
            in_group[perm[jp]] = -1; //flagged;
            out_group[perm[jp]] = new_group_idx;
            jp++;
        }
}

int saad_reordering(CSR& cmat, input_parameters &params, intT* out_group, int (*reorder_func)(CSR&, intT*, input_parameters&), bool (*sim_condition)(intT*, intT, intT*, intT, input_parameters&))
{

    intT* in_group = new intT[cmat.rows];
    reorder_func(cmat, in_group, params); //creates a first grouping for rows (hash-based)

    intT* perm = new intT[cmat.rows];
    sort_permutation(perm, in_group, cmat.rows);

    intT current_in_group = 0;
    intT current_out_group = 0;


    intT* group_structure; //holds the nz-structure of the current group 
    intT group_structure_nzcount;

    intT i, j;
    for (intT ip = 0; ip < cmat.rows; ip++)
    {
        i = perm[ip];
        if (in_group[i] != -1) assign_group(in_group, out_group, perm, cmat.rows, ip, current_out_group);

        intT last_checked = -2; //used to jump over already seen (but unassigned) groups;

        make_group_structure(group_structure, group_structure_nzcount, cmat.ja[i], cmat.nzcount[i], params);

        //check all (groups of) rows after i; 
        for (intT jp = ip + 1; jp < cmat.rows; jp++)
        {
            j = perm[jp];
            if (in_group[j] != -1 && in_group[j] != last_checked)
            {
                last_checked = in_group[j];

                if (sim_condition(group_structure, group_structure_nzcount, cmat.ja[j], cmat.nzcount[j], params))
                {
                    assign_group(in_group, out_group, perm, cmat.rows, jp, current_out_group);
                    update_group_structure(group_structure, group_structure_nzcount, cmat.ja[j], cmat.nzcount[j], params);
                }
            }
        }
        current_out_group++;
    }


    delete[] in_group;
    if (group_structure) delete[] group_structure;

    return 0;

}


int make_group_structure(intT* group_structure, intT &group_structure_nzcount, intT* cols_A, intT len_A, input_parameters &params)
{
    if (len_A == 0) return 0;

    if (params.reorder_algo == "saad")
    {
        group_structure = new intT[len_A];
        group_structure_nzcount = len_A;
        std::copy(cols_A, len_A, group_structure);
    }
    else if (params.reorder_algo == "saad_blocks")
    {
        intT block_size = params.algo_block_size;
        group_structure = new intT[len_A];
        
        intT current_block = -1;
        intT group_idx = 0;
        for (intT i = 0; i < len_A; i++)
        {
            current_block = cols_A[i] % block_size;
            group_structure[group_idx] = current_block;
            while (i < len_A && cols_A[i] % block_size == current_block) i++;
            group_idx++;
        }
        group_structure_nzcount = group_idx;
    }
    else
    {
        return 1;
    }
    return 0;
}

int update_group_structure(intT* group_structure, intT& group_structure_nzcount, intT* cols_A, intT len_A, input_parameters& params)
{

    intT block_size;
    if (params.reorder_algo == "saad") block_size = 1;
    else if (params.reorder_algo == "saad_blocks") block_size = params.algo_block_size;
    else return 1; //unknown reordering algorithm

    if (len_A + group_structure_nzcount == 0) return 0;

    intT* new_group_structure = new intT[len_A + group_structure_nzcount];
    intT group_idx = 0;
    intT new_group_idx = 0;
    intT j = 0;
    while (group_idx < group_structure_nzcount || j < len_A)
    {
        if (group_idx < group_structure_nzcount && group_structure[group_idx] < cols_A[j] / block_size)
        {
            new_group_structure[new_group_idx] = group_structure[group_idx];
            group_idx++;
            new_group_idx++;
        }
        else if (j < len_A && group_structure[group_idx] > cols_A[j] / block_size)
        {
            new_group_structure[new_group_idx] = cols_A[j] / block_size;
            new_group_idx++;
            intT current_block = cols_A[j] / block_size;
            while (j < len_A && cols_A[j] / block_size == current_block) j++;
        }
        else if (j < len_A && group_idx < group_structure_nzcount)
        {
            new_group_structure[new_group_idx] = group_structure[group_idx];
            new_group_idx++;
            group_idx++;
            while (j < len_A && cols_A[j] / block_size == current_block) j++;
        }
    }

    delete[] group_structure;
    group_structure = new_group_structure;
    group_structure_nzcount = new_group_idx;
    return 0;
}

int saad_reordering(CSR& cmat, input_parameters& params, intT* out_group)
{  
    if (params.reorder_algo == "saad")
        saad_reordering(cmat, params, out_group, hash_reordering, scalar_condition);
    else if (params.reorder_algo == "saad_blocks")
        saad_reordering(cmat, params, out_group, hash_reordering, scalar_block_condition);
    else
        std::cout << "UNKNONW ALGORITMH -->" << params.reorder_algo << "<-- in saad reordering" << std::endl;
}

bool scalar_condition(intT* cols_A, intT len_A, intT* cols_B, intT len_B, input_parameters &params)
{

    float eps = params.eps;
    if (len_A == 0 && len_B == 0) return true;
    if (len_A == 0 || len_B == 0) return false;


    intT count = 0;
    intT i = 0, j = 0;

    while (i < len_A && j < len_B)
    {
        if (cols_A[i] < cols_B[j]) i++;
        else if (cols_A[i] > cols_B[j]) j++;
        else
        {
            i++;
            j++;
            count++;
        }
    }

    if ((std::pow(count, 2) > std::pow(eps, 2) * len_A * len_B)) return true;
    else return false;

}

bool scalar_block_condition_csr(intT* cols_A, intT len_A, intT* cols_B, intT len_B, input_parameters &params)
{

    float eps = params.eps;
    intT block_size = params.algo_block_size;
    if (len_A == 0 && len_B == 0) return true;
    if (len_A == 0 || len_B == 0) return false;

    intT modA, modB;
    intT len_mod_A = 0;
    for(intT i = 0; i < len_A;)
    {
        len_mod_A++;
        modA = cols_A[i] / block_size;
        while (i < len_A && cols_A[i] / block_size == modA) i++;
    }

    intT len_mod_B = 0;
    for (intT i = 0; i < len_B;)
    {
        len_mod_B++;
        modB = cols_B[i] / block_size;
        while (i < len_B && cols_B[i] / block_size == modB) i++;
    }

    intT i = 0;
    intT j = 0;
    intT count = 0;
    while (i < len_A && j < len_B)
    {
        modA = cols_A[i] / block_size;
        modB = cols_B[j] / block_size;
        if (modA == modB)
        {
            count++;
            modA++;
            modB++;
        }
        while (i < len_A && cols_A[i]/block_size < modB) i++;
        while (j < len_B && cols_B[j]/block_size < modA) j++;
    }

    if ((std::pow(count, 2) > std::pow(eps, 2) * len_mod_A * len_mod_B)) return true;
    else return false;
}

bool scalar_block_condition(intT* group_structure, intT group_structure_nzcount, intT* cols_B, intT len_B, input_parameters& params)
{

    float eps = params.eps;
    intT block_size = params.algo_block_size;
    if (len_A == 0 && group_structure_nzcount == 0) return true;
    if (len_A == 0 || group_structure_nzcount == 0) return false;

    intT modB;
    intT len_mod_B = 0; //counts the size of cols_B when partitioned with block_size;
    for (intT j = 0; j < len_B;)
    {
        len_mod_B++;
        modB = cols_B[j] / block_size;
        while (j < len_B && cols_B[j] / block_size == modB) i++;
    }

    intT group_idx = 0;
    intT j = 0;
    intT count = 0;
    while (group_idx < group_structure_nzcount && j < len_B)
    {
        modA = group_structure[group_idx];
        modB = cols_B[j] / block_size;
        if (group_structure[group_idx] == modB)
        {
            count++;
            modA++;
            modB++;
        }

        while (group_idx < group_structure_nzcount && group_structure[group_idx] < modB) i++;
        while (j < len_B && cols_B[j] / block_size < modA) j++;
    }

    if ((std::pow(count, 2) > std::pow(eps, 2) * group_structure_nzcount * len_mod_B)) return true;
    else return false;
}


int group_to_VBS(CSR& cmat, intT* grouping, intT* compressed_dim_partition, intT nB, VBS& vbmat, int vbmat_blocks_fmt, int vbmat_entries_fmt)
{
    //create a VBS reordering the main (uncompressed) dimension of a CSR matrix according to a grouping

    intT rows = cmat.rows;
    intT cols = cmat.cols;
    intT main_dim = (cmat.fmt == 0) ? rows : cols;

    intT* perm = new intT[main_dim];

    sort_permutation(perm, grouping, main_dim); //find a permutation that sorts the main dimension according to the grouping

    intT grp_num;
    grp_num = count_groups(grouping, main_dim);

    //partition the main dimension
    intT* main_partition = new intT[grp_num + 1];
    grp_to_partition(grouping, main_dim, main_partition); 

    //create a permuted CSR
    CSR cmat_cpy;
    copy(cmat, cmat_cpy);

    permute_CSR(cmat_cpy, perm, cmat_cpy.fmt);

    intT* row_part;
    intT row_blocks;
    intT* col_part;
    intT col_blocks;

    if (cmat.fmt == 0)
    {
        row_part = main_partition;
        row_blocks = grp_num;
        col_part = compressed_dim_partition;
        col_blocks = nB;
    }
    else
    {
        col_part = main_partition;
        col_blocks = grp_num;
        row_part = compressed_dim_partition;
        row_blocks = nB;
    }

    convert_to_VBS(cmat_cpy,
        vbmat,
        row_blocks, row_part,
        col_blocks, col_part,
        vbmat_blocks_fmt, vbmat_entries_fmt);

    //cleaning

    cleanCSR(cmat_cpy);
    delete[] perm;
    delete[] main_partition;
}