# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import numpy as np
import matplotlib.pyplot as plt
import glob



def evaluate_blocking(grouping,nz_block_count,block_width):
    #get list of row_block_sizes from grouping;
    row_block_heights = []
    last_visited = -1
    for group in np.sort(np.array(grouping)):
        if group != last_visited:
            row_block_heights.append(1)
        row_block_heights[-1] += 1
        last_visited = group

    nztot = 0
    for count,block_height in zip(nz_block_count, row_block_heights):
        nztot += count*block_height*block_width

    return nztot, row_block_heights


def extract_data(filename):
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
                data = {var:(float(value) if value.isnumeric() else value) for var, value in zip(vars, line)}
            elif idx == 2:
                nzcount = [float(elem) for elem in line[1:]]
            elif idx == 3:
                grouping = [float(elem) for elem in line[1:]]
                complete = True
    if complete:
        return data,nzcount,grouping
    else:
        print(f"ERROR: MISSING DATA in {filename}")
        return {},[],[]



def check_constraints(data_dict, constraints_dict):
    for constraint in constraints_dict:
        if not constraint in data_dict:
            return False
        else:
            if constraints_dict[constraint] == data_dict[constraint]:
                return True
            else:
                return False 

def get_data_line(folder, constraints):
    datapoints = []
    for experiment in glob.glob(f"{folder}*/*.txt"):
        data,nzcount,grouping = extract_data(experiment)
        if check_constraints(data,constraints):
            nztot, row_block_heights = evaluate_blocking(grouping, nzcount, data["col_block_size"])
            data["nonzeros_padding"] = nztot
            data["avg_height"] = np.average(row_block_heights)
            datapoints.append(data)
            #print(data)
    return datapoints


def check_unique(l):
    return len(l) > len(set(l))

def plot_against(folder, x_name = "tau", y_name = "nonzeros_padding", constraints = {}):
    datapoints = get_data_line(folder, constraints)
    x_points = [data[x_name] for data in datapoints]
    y_points = [data[y_name] for data in datapoints]
    if not check_unique(x_points) or not check_unique(y_points):
        print("WARNING: plotting elements are not unique. Define better constraints")
    print(x_points, y_points)
    plt.plot(x_points, y_points)


constraints = {"blocking_algo": 0, 
               "col_block_size": 32,
                }

plot_against("results/ia", constraints = constraints)
