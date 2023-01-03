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
        std::string outfile_ = "";
        std::string exp_name_ = "";
        std::string reader_delimiter_ = " ";
        bool scramble_ = false;
        int sim_measure_ = 1;
        bool sim_use_groups_ = 1;
        bool sim_use_pattern_ = 1;
        bool pattern_only_ = 0;
        int block_repetitions_ = 1;
        float tau_ = 0.5;
        int block_size_ = 1;

        CLineReader(int argc, char* argv[])
        {
            ParseArgs(argc, argv);
        }

        bool ParseArgs(int argc, char* argv[])
        {
            char c_opt;
            while ((c_opt = getopt(argc, argv, "b:f:g:n:o:p:r:s:t:")) != -1)
            {
                switch(c_opt) 
                {
                    case 'b': block_size_ = std::stoi(optarg);              break;
                    case 'g': sim_use_groups_ = (std::stoi(optarg) == 1);   break;
                    case 'o': outfile_ = std::string(optarg);          break;
                    case 'p': sim_use_pattern_ = (std::stoi(optarg) == 1);  break;
                    case 'm': sim_measure_ = std::stoi(optarg);             break;
                    case 'n': exp_name_ = std::string(optarg);         break;
                    case 'f': filename_ = std::string(optarg);         break;
                    case 'r': block_repetitions_ = std::stoi(optarg);       break;
                    case 's': scramble_ = (std::stoi(optarg) == 1);         break;
                    case 't': tau_ = std::stof(optarg);                     break;

                }           
            }

            if ((filename_ == "")) 
            {
            std::cout << "No graph input specified." << std::endl;
            return false;
            }
        }
};