#pragma once
#include "matrices.h"

intT HammingDistance(intT* row_A, intT size_A, intT* row_B, intT size_B);
float JaccardDistance(intT* row_A, intT size_A, intT* row_B, intT size_B);
float CosineDistance(intT* row_A, intT size_A, intT* row_B, intT size_B);



intT HammingDistanceQuotient(intT* row_A, intT size_A, intT* row_B, intT size_B, intT stride);
intT JaccardDistanceQuotient(intT* row_A, intT size_A, intT* row_B, intT size_B, intT stride);

