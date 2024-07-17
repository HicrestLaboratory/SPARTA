#include "matrices.h"
#include "blocking.h"
#include "utilities.h"

#include <fstream>
#include <iostream>

using namespace std;

std::vector<intT> read_grouping_file(const std::string& filename) {
    std::vector<intT> grouping;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return grouping;
    }

    std::string line;
    while (std::getline(infile, line)) {
        try {
            intT number = std::stoi(line);
            grouping.push_back(number);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid number in file " << filename << " at line: " << line << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Error: Number out of range in file " << filename << " at line: " << line << std::endl;
        }
    }

    infile.close();
    return grouping;
}



int main(int argc, char* argv[])
{
    bool reorder;
    std::string input_file;
    std::string grouping_file;
    int block_size;
    bool symmetric_reorder = false;

    if (argc == 3) {
        reorder = false;
    }
    else if (argc == 5){
        reorder = true;
        grouping_file = argv[3];
        symmetric_reorder = std::stoi(argv[4]);
    }

    else{
        std::cerr << "Usage: " << argv[0] << "<input_matrix> <block_size> (OPTIONAL) <grouping file> <reorder-symmetric>" << std::endl;
        return 1;
    }

    input_file = argv[1];
    block_size = std::stoi(argv[2]);
    std::string delimiter = " ";
    bool pattern_only = false;
    MatrixFormat mat_fmt = mtx;
    bool symmetrize = false;

    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << input_file << std::endl;
        return 1;
    }

    CSR csr_matrix(infile);
    infile.close();

    if (reorder)
    {
        //READ GROUPING
        std::vector<intT> grouping = read_grouping_file(grouping_file);
        if (grouping.size() == csr_matrix.rows + 1) grouping.erase(grouping.begin());

        //PERMUTE MATRIX
        if (!symmetric_reorder){
            csr_matrix.reorder(grouping);
            //std::cerr << "Matrix successfully reordered (only rows) according to file: " << grouping_file << std::endl;
        }
        else 
        {
            csr_matrix.reorder2d(grouping);
            //std::cerr << "Matrix successfully reordered (2d) according to file: " << grouping_file << std::endl;
        }
    }

    BlockingEngine bEngine;
    bEngine.col_block_size = block_size;
    bEngine.row_block_size = block_size;
    bEngine.force_fixed_size = false;
    bEngine.grouping_result = FixedBlocking(csr_matrix, block_size);
    bEngine.CollectBlockingInfo(csr_matrix);
    std::cout << bEngine.VBR_nzcount << " " << bEngine.VBR_nzblocks_count << " " << bEngine.VBR_average_height << " " << bEngine.VBR_longest_row << std::endl;




    //    vector<intT> grouping = //LOAD GROUPING FROM FILE

}
