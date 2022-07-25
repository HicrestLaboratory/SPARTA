#include "similarity.h"
#include <vector>
#include <iostream>

using namespace std;

#pragma once
#include "matrices.h"

float HammingDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  if (size_A == 0 and size_B == 0) return 0;
  if (size_A == 0 or size_B == 0) return max(size_A,size_B);
  intT i = 0;
  intT j = 0;

  cout << "NOT IMPLEMENTED YET" << endl;
  return 0;
}

float JaccardDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  cout << "NOT IMPLEMENTED YET" << endl;
  return 0;
}

float CosineDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  cout << "NOT IMPLEMENTED YET" << endl;
  return 0;
}
