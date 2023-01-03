#include "blocking.h"
#include "matrices.h"
#include "utilities.h"
#include "input.h"
#include <math.h> //ceil, floor
#include <vector>
#include <iostream>

using namespace std;
std::vector<intT> IterativeBlockingPattern(const CSR& cmat, float tau, distFuncGroup distanceFunction,intT block_size, bool use_size, bool use_pattern, intT &comparison_counter, intT &merge_counter, float &timer)
{
    vector<intT> grouping(cmat.rows, -1); //flag each rows as ungrouped (-1)

    //main loop. Takes an ungrouped row (i) as a seed for a new block.
    for (intT i = 0; i < cmat.rows; i++)
    {
        if (grouping[i] == -1)
        {
            vector<intT> pattern;
            intT current_group_size = 1;
            grouping[i] = i; // the group is named the same as the seed row.
            pattern.insert(pattern.end(), &cmat.ja[i][0], &cmat.ja[i][cmat.nzcount[i]]); //Initialize the pattern with the seed entries.

            //inner loop, compare each subsequent row with the current pattern
            for (intT j = i + 1; j < cmat.rows; j++)
            {
                if (grouping[j] == -1)
                {
                    comparison_counter++;
                    float dist = distanceFunction(pattern, current_group_size, cmat.ja[j], cmat.nzcount[j], 1, block_size);
                    if (dist < tau)
                    {
                        merge_counter++;
                        grouping[j] = i;
                        if (use_pattern)
                            pattern = merge_rows(pattern, cmat.ja[j], cmat.nzcount[j]); //update the pattern through union with the new row
                        if (use_size)
                            current_group_size++;
                    }
                }
            }
        }
    }
    return grouping;
}

std::vector<intT> BlockingEngine::ObtainPartition(const CSR& cmat)
{
    //run the blocking function and store statistics
    comparison_counter = 0;
    merge_counter = 0;
    timer = 0;
    return IterativeBlockingPattern(cmat, tau, comparator, block_size, use_groups, use_pattern, comparison_counter, merge_counter, timer);
}

BlockingEngine::BlockingEngine(CLineReader &cline)
{       
  tau = cline.tau_;
  block_size = cline.block_size_;
  use_groups = cline.sim_use_groups_;
  use_pattern = cline.sim_use_pattern_;
  SetComparator(cline.sim_measure_);
}

void BlockingEngine::SetComparator(int choice)
{
    switch(choice)
    {
        case 0:
            comparator = HammingDistanceGroup;
            break;
        case 1: 
            comparator = JaccardDistanceGroup;
            break;
    }
}

float HammingDistanceGroup(vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size)
{
  bool count_zeros = 0;
  intT size_A = row_A.size();
  if (size_A == 0 and size_B == 0) return 0;

  intT blocksize_A = ceil(size_A/block_size);
  intT blocksize_B = ceil(size_B/block_size);
  
  if (size_A == 0 or size_B == 0) return max(blocksize_A*group_size_A,blocksize_B*group_size_B);

  intT add_to_count_A, add_to_count_B;
  if (count_zeros)
  {
    add_to_count_A = group_size_B;
    add_to_count_B = group_size_A;
  }
  else
  {
    add_to_count_A = group_size_A;
    add_to_count_B = group_size_B;
  }

  auto blockpos_A = [&](intT i){return floor(row_A[i]/block_size);}; 
  auto blockpos_B = [&](intT j){return floor(row_B[j]/block_size);};


  intT i = 0;
  intT j = 0;
  intT count = 0;
  intT pos_A;
  intT pos_B;
  
  while (i < size_A && j < size_B)
  {
    pos_A = blockpos_A(i);
    pos_B = blockpos_B(j);

    if (pos_A < pos_B)
    {
      count += add_to_count_A;
      while(i < size_A && blockpos_A(i) == pos_A) i++;
    }
    else if (pos_A > pos_B)
    {
      count += add_to_count_B;
      while(j < size_B && blockpos_B(j) == pos_B) j++;
    }
    else
    {
      while(i < size_A && blockpos_A(i) == pos_A) i++;
      while(j < size_B && blockpos_B(j) == pos_B) j++;
    }
  }

  while (i < size_A)
  { 
    pos_A = blockpos_A(i);
    count += add_to_count_A;
    while(i < size_A && blockpos_A(i) == pos_A) i++;
  }

  while (j < size_B)
  { 
    pos_B = blockpos_B(j);
    count += add_to_count_B;
    while(j < size_B && blockpos_B(j) == pos_B) j++;
  }

  return count;
}


float JaccardDistanceGroup(vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size)
{
  intT size_A = row_A.size();
  if (size_A == 0 && size_B == 0) return 0;
  if (size_A == 0 || size_B == 0) return 1;

  float h = HammingDistanceGroup(row_A, group_size_A, row_B, size_B, group_size_B, block_size);
  return 2*h/(size_A*group_size_A + size_B*group_size_B + h);
}