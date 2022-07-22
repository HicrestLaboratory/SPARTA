#include "similarity.h"
#include <vector>

using namespace std;

#pragma once
#include "matrices.h"

float HammingDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  if (size_A == 0 and size_B == 0) return 0;
  if (size_A == 0 or size_B == 0) return max(size_A,size_B);
  intT i = 0;
  intT j = 0;
  While (
  
}

float JaccardDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
 
}

float CosineDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  
}
