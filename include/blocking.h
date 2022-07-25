#pragma once
#include <vector>
#include <string>
#include "matrices.h"
#include "similarity.h"

typedef int (*distFunc)(intT*,intT,intT*,intT);

void IterativeBlockingJaccard(const CSR& cmat, intT* grouping);
void IterativeBlockingGeneral(const CSR& cmat, intT* grouping, float tau = 0.5, distFunc distanceFunction = &JaccardDistance);