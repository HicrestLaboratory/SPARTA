#pragma once

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unistd.h> //getopt, optarg

class CLineReader
{
    public:
        std::string filename_ = "data/TEST_matrix_weighted.txt";
        std::string outfile_ = "results/TEST_results.txt";
        std::string exp_name_ = "";
        std::string reader_delimiter_ = " ";
        
        bool scramble_ = false;
        bool sim_use_groups_ = 1;
        bool sim_use_pattern_ = 1;
        bool pattern_only_ = 0;
        
        int sim_measure_ = 1;
        int block_repetitions_ = 1;
        int block_size_ = 1;
        int verbose_ = 0;

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
            std::cout << "scramble_: " <<  scramble_ << std::endl;
            std::cout << "sim_measure_: " <<  sim_measure_ << std::endl;
            std::cout << "sim_use_groups_: " <<  sim_use_groups_ << std::endl;
            std::cout << "sim_use_pattern_: " <<  sim_use_pattern_ << std::endl;
            std::cout << "pattern_only_: " <<  pattern_only_ << std::endl;
            std::cout << "block_repetitions_: " <<  block_repetitions_ << std::endl;
            std::cout << "tau_: " <<  tau_ << std::endl;
            std::cout << "block_size_: " <<  block_size_ << std::endl;
            std::cout << "verbose_: " <<  verbose_ << std::endl;
            std::cout << "___________________" << std::endl;; 
        }


        void ParseArgs(int argc, char* argv[])
        {
            char c_opt;
            while ((c_opt = getopt(argc, argv, "b:f:g:n:o:p:P:r:s:t:v:")) != -1)
            {
                switch(c_opt) 
                {
                    case 'b': block_size_ = std::stoi(optarg);              break;
                    case 'g': sim_use_groups_ = (std::stoi(optarg) == 1);   break;
                    case 'o': outfile_ = std::string(optarg);          break;
                    case 'p': sim_use_pattern_ = (std::stoi(optarg) == 1);  break;
                    case 'P': pattern_only_ = (std::stoi(optarg) == 1);  break;
                    case 'm': sim_measure_ = std::stoi(optarg);             break;
                    case 'n': exp_name_ = std::string(optarg);         break;
                    case 'f': filename_ = std::string(optarg);         break;
                    case 'r': block_repetitions_ = std::stoi(optarg);       break;
                    case 's': scramble_ = (std::stoi(optarg) == 1);         break;
                    case 't': tau_ = std::stof(optarg);                     break;
                    case 'v': verbose_ = std::stoi(optarg);                     break;
                }           
            }

            if ((filename_ == "")) 
            {
            std::cout << "No graph input specified." << std::endl;
            }
        }
};