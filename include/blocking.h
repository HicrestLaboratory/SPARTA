#pragma once
#include <vector>
#include <string>
#include "matrices.h"
#include "similarity.h"

typedef float (*distFunc)(intT*,intT,intT*,intT);
typedef float (*distFuncGroup)(std::vector<intT>,intT,intT*,intT,intT,intT);

std::vector<intT> IterativeBlockingGeneral(const CSR& cmat, float tau = 0.5, distFunc distanceFunction = &JaccardDistance);
std::vector<intT> IterativeBlockingPattern(const CSR& cmat, float tau, distFuncGroup distanceFunction, intT block_size);
