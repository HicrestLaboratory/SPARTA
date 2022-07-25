#include "similarity.h"
#include <vector>
#include <iostream>
#include <math.h> //ceil

using namespace std;

#pragma once
#include "matrices.h"


intT HammingDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  if (size_A == 0 and size_B == 0) return 0;
  if (size_A == 0 or size_B == 0) return max(size_A,size_B);

  intT i = 0;
  intT j = 0;
  intT count = 0;

  while (i < size_A && j < size_B)
  {
    if (row_A[i] < row_B[j])
    {
      count++;
      i++;
    }
    else if (row_A[i] > row_B[j])
    {
      count++;
      j++;
    }
    else
    {
      i++;
      j++;
    }
  }
  if (i < size_A) count += size_A - i;
  if (j < size_B) count += size_B - j;

  return count;
}

float JaccardDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  if (size_A == 0 && size_B == 0) return 0;
  if (size_A == 0 || size_B == 0) return 1;

  float h = HammingDistance(row_A, size_A, row_B, size_B);
  return 2*h/(size_A + size_B + h);
}

float CosineDistance(intT* row_A, intT size_A, intT* row_B, intT size_B)
{
  cout << "NOT IMPLEMENTED YET" << endl;
  return 0;
}


intT HammingDistanceQuotient(intT* row_A, intT size_A, intT* row_B, intT size_B, intT stride)
{
  if (size_A == 0 and size_B == 0) return 0;

  intT blocksize_A = ceil(size_A/stride);
  intT blocksize_B = ceil(size_B/stride);

  if (size_A == 0 or size_B == 0) return max(blocksize_A,blocksize_B);

  intT i = 0;
  intT j = 0;
  intT count = 0;

  auto blockpos_A = [&](){return floor(row_A[i]/stride)}; 
  auto blockpos_B = [&](){return floor(row_B[j]/stride)};

  while (i < size_A || j < size_B)
  {
    intT pos_A = blockpos_A();
    intT pos_B = blockpos_B();

    if (pos_A < pos_B || j >= size_B)
    {
      count++;
      while(i < size_A && blockpos_A() == pos_A) i++;
    }
    else if (pos_A > pos_B || i >= size_A)
    {
      count++;
      while(j < size_B && blockpos_B() == pos_B) j++;
    }
    else
    {
      while(i < size_A && blockpos_A() == pos_A) i++;
      while(j < size_B && blockpos_B() == pos_B) j++;
    }
  }

  return count;
}
