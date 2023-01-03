#pragma once
#include <vector>
#include <string>
#include "matrices.h"
#include "input.h"


typedef float (*distFunc)(intT*,intT,intT*,intT);
typedef float (*distFuncQuot)(intT*,intT,intT*,intT, intT );
typedef float (*distFuncGroup)(std::vector<intT>,intT,intT*,intT,intT,intT);

class BlockingEngine
{
    public:
        float tau = 0.5;
        int block_size = 1;
        bool use_groups = 0;
        bool use_pattern = 1;
        float timer = 0;
        intT comparison_counter = 0;
        intT merge_counter = 0;
        distFuncGroup comparator;
        std::vector<intT> ObtainPartition(const CSR& cmat);
        void SetComparator(int choice);

        BlockingEngine(){};
        BlockingEngine(CLineReader &cline);
};

std::vector<intT> IterativeBlockingPattern(const CSR& cmat, float tau, distFuncGroup distanceFunction,intT block_size, bool use_size, bool use_pattern,  intT &comparison_counter, intT &merge_counter, float &timer);

float HammingDistanceGroup(std::vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size);
float JaccardDistanceGroup(std::vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size);

