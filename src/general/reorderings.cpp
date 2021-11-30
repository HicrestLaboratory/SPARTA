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
    if (len_A == 0) return true;
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

    intT block_size = params.reorder_algo == "saad" ? 1 : params.algo_block_size;

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
            if (!check_same_pattern(cmat.ja[curr], cmat.nzcount[curr], cmat.ja[prev], cmat.nzcount[prev], block_size, 0)) current_group++;
        }
        groups[curr] = current_group;
    }
}

intT assign_group(intT* in_group, intT* out_group, intT* perm, intT len, intT jp, intT new_group_idx)
{
    intT group_size = 0;
    intT current_in_group = in_group[perm[jp]];
    while (jp < len && in_group[perm[jp]] == current_in_group)
        {
            group_size++,
            in_group[perm[jp]] = -1; //flagged;
            out_group[perm[jp]] = new_group_idx;
            jp++;
        }
    return group_size;
}

int saad_reordering(CSR& cmat, input_parameters &params, intT* out_group, int (*reorder_func)(CSR&, intT*, input_parameters&), bool (*sim_condition)(group_structure& group_struct, intT*, intT, intT, input_parameters&), reorder_info& info)
{

    intT* in_group = new intT[cmat.rows];
    reorder_func(cmat, in_group, params); //creates a first grouping for rows (hash-based)

    intT* perm = new intT[cmat.rows];
    sort_permutation(perm, in_group, cmat.rows);

    intT current_out_group = 0;
   
    intT second_group_size;

    intT i, j;

    intT j_tmp;

    for (intT ip = 0; ip < cmat.rows; ip++)
    {

        i = perm[ip];

        group_structure group_struct; //holds the nz-structure of the current group 

        if (in_group[i] != -1)
        {

            assign_group(in_group, out_group, perm, cmat.rows, ip, current_out_group);

            intT last_checked = -2; //used to jump over already seen (but unassigned) groups;

            make_group_structure(group_struct, cmat.ja[i], cmat.nzcount[i], params);

            //check all (groups of) rows after i; 
            for (intT jp = ip + 1; jp < cmat.rows; jp++)
            {
                j = perm[jp];
                if (in_group[j] != -1 && in_group[j] != last_checked)
                {
                    last_checked = in_group[j];

                    //count size of the group to be compared;
                    j_tmp = jp+1;
                    second_group_size = 1;
                    while (j_tmp < cmat.rows && in_group[perm[j_tmp]] == last_checked)
                    {
                        second_group_size++;
                        j_tmp++;
                    }
                    //---


                    info.comparisons++;
                    if (sim_condition(group_struct, cmat.ja[j], cmat.nzcount[j], second_group_size, params))
                    {
                        assign_group(in_group, out_group, perm, cmat.rows, jp, current_out_group);
                        update_group_structure(group_struct, cmat.ja[j], cmat.nzcount[j], second_group_size, params);
                    }
                }
            }
            current_out_group++;
        }

        info.skipped += group_struct.skipped;
        group_struct.clean();
    }


    delete[] in_group;

    return 0;

}

int saad_reordering(CSR& cmat, input_parameters& params, intT* out_group, reorder_info& info)
{  
    if (params.reorder_algo == "saad")
        saad_reordering(cmat, params, out_group, hash_reordering, scalar_condition, info);
    else if (params.reorder_algo == "saad_blocks")
        saad_reordering(cmat, params, out_group, hash_reordering, scalar_block_condition, info);
    else
        std::cout << "UNKNONW ALGORITMH -->" << params.reorder_algo << "<-- in saad reordering" << std::endl;
}

int saad_reordering(CSR& input_cmat, VBS& output_vbmat, intT algo_block_size, int vbmat_blocks_fmt, int vbmat_entries_fmt, input_parameters& params, reorder_info& info)
{
    vbmat_blocks_fmt = 1;
    vbmat_entries_fmt = 1;
    intT block_cols = std::ceil((float)input_cmat.cols / algo_block_size);

    //prepare the column partition
    intT* col_part = new intT[block_cols + 1];
    partition(col_part, 0, input_cmat.cols, algo_block_size);

    //run the reordering algo
    intT* hash_groups = new intT[input_cmat.rows];
    saad_reordering(input_cmat, params, hash_groups, info);

    //create the block matrix
    group_to_VBS(input_cmat, hash_groups, col_part, block_cols, output_vbmat, vbmat_blocks_fmt, vbmat_entries_fmt);

    delete[] hash_groups;

}


bool scalar_condition(group_structure& group_struct, intT* cols_B, intT len_B, intT group_size_B, input_parameters& params)
{
    float eps = params.eps;
    if (group_struct.len == 0 && len_B == 0) return true;
    if (group_struct.len == 0 || len_B == 0) return false;

    intT count = 0;
    intT i = 0, j = 0;

    while (i < group_struct.len && j < len_B)
    {
        if (group_struct.structure[i] < cols_B[j]) i++;
        else if (group_struct.structure[i] > cols_B[j]) j++;
        else
        {
            i++;
            j++;
            count++;
        }
    }

    bool result;
    if (params.similarity_func == "hamming") result = group_struct.len + len_B - (2 * count) < eps * params.A_cols;
    else if (params.similarity_func == "scalar") result = (std::pow(count, 2) > std::pow(eps, 2) * group_struct.len * len_B);
    else if (params.similarity_func == "jaccard") result = (1.0 * count) / (group_struct.len + len_B - count) > eps;

    return result;

}

///DEPRECATED
/*bool scalar_block_condition(intT* group_structure, intT group_structure_nzcount, intT group_size, intT* cols_B, intT len_B, intT group_size_B, input_parameters& params)
{

    float eps = params.eps;
    intT block_size = params.algo_block_size;
    if (len_B == 0 && group_structure_nzcount == 0) return true;
    if (len_B == 0 || group_structure_nzcount == 0) return false;

    intT modB;
    intT len_mod_B = 0; //counts the size of cols_B when partitioned with block_size;
    for (intT j = 0; j < len_B;)
    {
        len_mod_B++;
        modB = cols_B[j] / block_size;
        while (j < len_B && cols_B[j] / block_size == modB) j++;
    }

    intT group_idx = 0;
    intT j = 0;
    intT count = 0;
    intT modA;
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

        while (group_idx < group_structure_nzcount && group_structure[group_idx] < modB) group_idx++;
        while (j < len_B && cols_B[j] / block_size < modA) j++;
    }


    bool result;
    if (params.similarity_func == "hamming") result = len_mod_B + group_structure_nzcount - (2 * count) < eps * params.A_cols;
    else if (params.similarity_func == "scalar") result = (std::pow(count, 2) > std::pow(eps, 2) * group_structure_nzcount * len_mod_B);
    else if (params.similarity_func == "jaccard") result = (1.0 * count) / (len_mod_B + group_structure_nzcount - count) > eps;
    else if (params.similarity_func == "density") result = ((float)(group_size * group_structure_nzcount + group_size_B * len_B) / ((group_structure_nzcount + len_B - count) * (group_size + group_size_B))) > eps;

    return result;
}
*/

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

bool scalar_block_condition(group_structure& group_struct, intT* cols_B, intT len_B, intT group_size_B, input_parameters& params)
{
    float eps = params.eps;
    intT block_size = params.algo_block_size;
    if (len_B == 0 && group_struct.len == 0) return true;
    if (len_B == 0 || group_struct.len == 0) return false;

    intT modB;
    intT len_mod_B = 0; //counts the size of cols_B when partitioned with block_size;
    for (intT j = 0; j < len_B;)
    {
        len_mod_B++;
        modB = cols_B[j] / block_size;
        while (j < len_B && cols_B[j] / block_size == modB) j++;
    }

    intT group_idx = 0;
    intT j = 0;
    intT count = 0;
    intT modA;
    while (group_idx < group_struct.len && j < len_B)
    {
        modA = group_struct.structure[group_idx];
        modB = cols_B[j] / block_size;
        if (group_struct.structure[group_idx] == modB)
        {
            count++;
            modA++;
            modB++;
        }

        while (group_idx < group_struct.len && group_struct.structure[group_idx] < modB) group_idx++;
        while (j < len_B && cols_B[j] / block_size < modA) j++;
    }

    bool result;
    if (params.similarity_func == "hamming") result = len_mod_B + group_struct.len- (2 * count) < eps * params.A_cols;
    else if (params.similarity_func == "scalar") result = (std::pow(count, 2) > std::pow(eps, 2) * group_struct.len * len_mod_B);
    else if (params.similarity_func == "jaccard") result = (1.0 * count) / (len_mod_B + group_struct.len - count) > eps;

    if (result && params.merge_limit != 0)
    {
        float limit_factor;
        if (params.merge_limit == -1)
        {
            limit_factor = 1 + 2. * ((1. - eps) / (3. - eps)); // the limit on the column number relative increase;
        }
        else
        {
            limit_factor = 1 + params.merge_limit;
        }

        if ((group_struct.len + len_mod_B - count) > limit_factor * group_struct.original_columns) //checks that the new number of columns is smaller than the bound
        {
            result = false;
            group_struct.skipped++; //count the number of lines skipped because of the merge bound
        }

    }
    return result;
}
int update_group_structure(group_structure& group_struct, intT* cols_A, intT len_A, intT A_group_size, input_parameters& params)
{
    //merges a blocked compressed row (group_structure) and a compressed row (cols_A) to update the blocked compressed row. 

    group_struct.group_size += A_group_size;

    if (params.hierarchic_merge == 0) return 0;

    intT block_size;
    if (params.reorder_algo == "saad") return 0; //do not use advanced merging for saad algorithm
    else if (params.reorder_algo == "saad_blocks") block_size = params.algo_block_size;
    else return 1; //unknown reordering algorithm

    if (len_A + group_struct.len == 0) return 0;


    intT* new_group_structure = new intT[len_A + group_struct.len];
    intT group_idx = 0;
    intT new_group_idx = 0;
    intT j = 0;
    intT current_block = 0;
    while (group_idx < group_struct.len && j < len_A)
    {
        if (group_struct.structure[group_idx] < cols_A[j] / block_size)
        {
            new_group_structure[new_group_idx] = group_struct.structure[group_idx];
            group_idx++;
            new_group_idx++;
        }
        else if (group_struct.structure[group_idx] > cols_A[j] / block_size)
        {
            new_group_structure[new_group_idx] = cols_A[j] / block_size;
            new_group_idx++;
            current_block = cols_A[j] / block_size;
            while (j < len_A && cols_A[j] / block_size == current_block) j++;
        }
        else if (group_idx < group_struct.len)
        {
            new_group_structure[new_group_idx] = group_struct.structure[group_idx];
            new_group_idx++;
            group_idx++;
            current_block = cols_A[j] / block_size;
            while (j < len_A && cols_A[j] / block_size == current_block) j++;
        }
    }

    while (group_idx < group_struct.len)
    {
        new_group_structure[new_group_idx] = group_struct.structure[group_idx];
        new_group_idx++;
        group_idx++;
    }

    while (j < len_A)
    {
        current_block = cols_A[j] / block_size;
        new_group_structure[new_group_idx] = current_block;
        new_group_idx++;
        while (j < len_A && cols_A[j] / block_size == current_block) j++;
    }

    if (group_struct.structure) delete[] group_struct.structure;
    group_struct.structure = new_group_structure;
    group_struct.len = new_group_idx;

    return 0;
}
int make_group_structure(group_structure& group_struct, intT* cols_A, intT len_A, input_parameters& params)
{
    if (len_A == 0)
    {
        group_struct.len = 0;
        return 0;
    }
    if (params.reorder_algo == "saad")
    {
        group_struct.structure = new intT[len_A];
        group_struct.len = len_A;
        std::copy(cols_A, cols_A + len_A, group_struct.structure);
    }
    else if (params.reorder_algo == "saad_blocks")
    {

        intT block_size = params.algo_block_size;
        group_struct.structure = new intT[len_A];

        intT current_block;
        intT group_idx = 0;
        for (intT i = 0; i < len_A; i++)
        {
            current_block = cols_A[i] / block_size;
            group_struct.structure[group_idx] = current_block;
            while (i < len_A && cols_A[i] / block_size == current_block) i++;
            group_idx++;
        }
        group_struct.len = group_idx;
        group_struct.original_columns = group_idx;
    }
    else
    {
        return 1;
    }
    return 0;
}
