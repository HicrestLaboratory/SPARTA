#pragma once
#include "matrices.h"

intT HammingDistance(intT* row_A, intT size_A, intT* row_B, intT size_B);
float JaccardDistance(intT* row_A, intT size_A, intT* row_B, intT size_B);

intT HammingDistanceQuotient(intT* row_A, intT size_A, intT* row_B, intT size_B, intT block_size);
float JaccardDistanceQuotient(intT* row_A, intT size_A, intT* row_B, intT size_B, intT block_size);

intT HammingDistanceGroup(intT* row_A, intT size_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size, bool count_zeros = true);
float JaccardDistanceGroup(intT* row_A, intT size_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size);
