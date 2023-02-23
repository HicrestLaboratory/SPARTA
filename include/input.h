#pragma once

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       //time
#include <random>       // random
#include <iostream>
#include <fstream>
#include <unistd.h> //getopt, optarg

class CLineReader
{
    public:
        std::string filename_ = "data/TEST_matrix_weighted.el";
        std::string outfile_ = "results/TEST_results.txt";
        std::string exp_name_ = "";
        std::string reader_delimiter_ = " ";
        
        bool sim_use_groups_ = 0;
        bool sim_use_pattern_ = 1;
        bool pattern_only_ = 0;
        bool force_fixed_size = 0;

        int blocking_algo_ = 3; // 0 for iterative; 1 for structured; 2 for fixed, 3 for iterative_clocked
        int seed_ = 0;
        int sim_measure_ = 1;
        int reorder_ = 0; //-1 for ascending. 0 for nothing. 1 for descending. 2 for scramble
        int col_block_size_ = 1;
        int row_block_size_ = 1;

        int verbose_ = 1;
        int warmup_ = 1; //how many warmup multiplications
        int exp_repetitions_ = 1;

        float tau_ = 0.5;

        CLineReader(int argc, char* argv[])
        {
            ParseArgs(argc, argv);
        }

        void print()
        {
            std::cout << "INPUT PARAMETERS:" << std::endl;
            std::cout << "filename: " << filename_ << std::endl;
            std::cout << "outfile: " <<  outfile_ << std::endl;
            std::cout << "exp_name: " <<  exp_name_ << std::endl;
            std::cout << "reader_delimiter_: " <<  reader_delimiter_ << std::endl;
            std::cout << "sim_measure_: " <<  sim_measure_ << std::endl;
            std::cout << "blocking_algo_: " <<  blocking_algo_ << std::endl;
            std::cout << "force_fixed_size: " <<  force_fixed_size << std::endl;
            std::cout << "sim_use_groups_: " <<  sim_use_groups_ << std::endl;
            std::cout << "sim_use_pattern_: " <<  sim_use_pattern_ << std::endl;
            std::cout << "pattern_only_: " <<  pattern_only_ << std::endl;
            std::cout << "reorder_by_degree_: " <<  reorder_ << std::endl;
            std::cout << "tau_: " <<  tau_ << std::endl;
            std::cout << "col_block_size_: " <<  col_block_size_ << std::endl;
            std::cout << "row_block_size_: " <<  row_block_size_ << std::endl;
            std::cout << "verbose_: " <<  verbose_ << std::endl;
            std::cout << "seed_: " <<  seed_ << std::endl; //-1 for random
            std::cout << "warmup_: " <<  warmup_ << std::endl;
            std::cout << "exp_repetitions_: " <<  exp_repetitions_ << std::endl;

            std::cout << "___________________" << std::endl;; 
        }


        void ParseArgs(int argc, char* argv[])
        {
            char c_opt;
            while ((c_opt = getopt(argc, argv, "a:b:B:f:F:g:m:n:o:p:P:r:s:t:v:w:x:")) != -1)
            {
                switch(c_opt) 
                {
                    case 'a': blocking_algo_ = std::stoi(optarg);                   break;
                    case 'b': col_block_size_ = std::stoi(optarg);                      break;
                    case 'B': row_block_size_ = std::stoi(optarg);                      break;
                    case 'f': filename_ = std::string(optarg);                      break;
                    case 'F': force_fixed_size = (std::stoi(optarg) == 1);                      break;
                    case 'g': sim_use_groups_ = (std::stoi(optarg) == 1);           break;
                    case 'o': outfile_ = std::string(optarg);                       break;
                    case 'p': sim_use_pattern_ = (std::stoi(optarg) == 1);          break;
                    case 'P': pattern_only_ = (std::stoi(optarg) == 1);             break;
                    case 'm': sim_measure_ = std::stoi(optarg);                     break;
                    case 'n': exp_name_ = std::string(optarg);                      break;
                    case 'r': reorder_ = std::stoi(optarg);                         break;
                    case 's': seed_ = std::stoi(optarg);                            break;
                    case 't': tau_ = std::stof(optarg);                             break;
                    case 'v': verbose_ = std::stoi(optarg);                         break;
                    case 'w': warmup_ = std::stoi(optarg);                          break;
                    case 'x': exp_repetitions_ = std::stoi(optarg);                 break;

                }           
            }

            if (seed_ != 0)
            {
                std::srand(seed_);
            }
            else
            {
                //TODO better randomness if needed
                std::random_device rd;
                std::srand(rd());
            }

            if ((filename_ == "")) 
            {
            std::cout << "No graph input specified." << std::endl;
            }
        }
};