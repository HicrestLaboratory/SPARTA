# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import numpy as np
import matplotlib.pyplot as plt


def evaluate_blocking(grouping,nz_block_count,block_width):
    #get list of row_block_sizes from grouping;
    row_block_heights = np.array()
    current = -1
    for group in np.sort(np.array(grouping)):
        if group != current:
            row_block_heights.append(1)
        row_block_heights[-1] += 1

    nztot = 0
    for count,block_height in zip(nz_block_count, row_block_heights):
        nztot += count*block_height*block_width

    return nztot, row_block_heights


def evaluate_blocking(filename):
    #example: 
    #matrix,rows,cols,nonzeros,blocking_algo,tau,row_block_size,col_block_size,use_pattern,sim_use_groups,sim_measure,exp_name,time_to_block,merge_counter,comparison_counter,
    #data/TEST_matrix_weighted.el,9,9,13,0,0.500000,1,1,1,1,1,,47.000000,2,31,
    #NZCOUNT,1,3,3,4,0,1,1,
    #GROUPING,0,1,2,3,4,5,4,4,8,

    with open(filename) as infile:
        complete = False
        for idx, line in enumerate(infile):
            line = line.split(",")[:-1]
            if idx == 0: 
                vars = line
            elif idx == 1:
                if len(line) != len(vars):
                    print(f"ERROR: CSV ENTRIES MISMATCH IN {filename}")
                    break
                data = {var:value for var, value in zip(vars, line)}
            elif idx == 2:
                nzcount = line[1:]
            elif idx == 3:
                grouping = line[1:]
                complete = True
    if complete:
        return data,nzcount,grouping
    else:
        print(f"ERROR: MISSING DATA in {filename}")
        return {},[],[]
        

def make_images(folder):


print(evaluate_blocking("results/TEST_results.txt"))