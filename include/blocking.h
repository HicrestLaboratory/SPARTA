#pragma once
#include <vector>
#include <string>
#include "matrices.h"
#include "input.h"
#include <fstream>


typedef float (*distFunc)(intT*,intT,intT*,intT);
typedef float (*distFuncQuot)(intT*,intT,intT*,intT, intT );
typedef float (*distFuncGroup)(std::vector<intT>,intT,intT*,intT,intT,intT);

enum BlockingType {iterative, iterative_pattern, fixed_size};


class BlockingEngine
{
    public:
        float tau = 0.5;

        intT col_block_size = 1;
        intT row_block_size = 1; //only used for fixed size blocking

        float timer = 0;
        
        bool use_groups = false;
        bool use_pattern = true;
        
        int structured_m = 2;
        int structured_n = 4;

        BlockingType blocking_algo = iterative;

        //bool blocking_completed = false;
        intT comparison_counter = 0;
        intT merge_counter = 0;
        distFuncGroup comparator;
        std::vector<intT> grouping_result;
        intT VBR_nzcount = 0;
        intT VBR_nzblocks_count = 0;
        float VBR_average_height = 0;

        std::vector<intT> GetGrouping(const CSR& cmat);
        void CollectBlockingInfo(const CSR& cmat);
        
        void SetComparator(int choice);
        
        void print();

        BlockingEngine(){};
        BlockingEngine(CLineReader &cline);
};

std::vector<intT> IterativeBlockingPattern(const CSR& cmat, float tau, distFuncGroup distanceFunction,intT block_size, bool use_size, bool use_pattern,  intT &comparison_counter, intT &merge_counter, float &timer);
std::vector<intT> IterativeBlockingPatternMN(const CSR& cmat, float tau, distFuncGroup distanceFunction,intT block_size, bool use_size, bool use_pattern, int structured_m, int structured_n, intT &comparison_counter, intT &merge_counter, float &timer);
std::vector<intT> FixedBlocking(const CSR& cmat, intT row_block_size);


float HammingDistanceGroup(std::vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size);
float JaccardDistanceGroup(std::vector<intT> row_A, intT group_size_A, intT* row_B, intT size_B, intT group_size_B, intT block_size);

