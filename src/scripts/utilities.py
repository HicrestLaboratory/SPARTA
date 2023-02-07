# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import operator



def evaluate_blocking(grouping,nz_block_count,block_width):
    #get list of row_block_sizes from grouping;
    row_block_heights = []
    last_visited = -1
    for group in np.sort(np.array(grouping)):
        if group != last_visited:
            row_block_heights.append(0)
        row_block_heights[-1] += 1
        last_visited = group

    nztot = 0
    for count,block_height in zip(nz_block_count, row_block_heights):
        nztot += count*block_height*block_width

    return nztot, row_block_heights


def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
    
def extract_data(filename):
    #example: 
    #matrix,rows,cols,nonzeros,blocking_algo,tau,row_block_size,col_block_size,use_pattern,sim_use_groups,sim_measure,reorder,exp_name,time_to_block,merge_counter,comparison_counter,
    #data/real_world//ia-wikiquote-user-edits-nodup.el,21608,94757,238714,0,0.010000,1,16,1,0,1,0,blocking_G_ia-wikiquote-user-edits-nodup_b_16_a_0_t_0.01_p_1_g_0_r_0,243922960.000000,10345,126953437,
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
                data = {var:(float(value) if isfloat(value) else value) for var, value in zip(vars, line)}
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
            print(f"WARNING: checking non existent property: {constraint} ")
            return False
        else:
            if constraints_dict[constraint] != data_dict[constraint]:
                return False 
    return True


def get_data_line(folder, constraints):
    datapoints = []
    for experiment in glob.glob(f"{folder}*/*.txt"):
        data,nzcount,grouping = extract_data(experiment)
        if check_constraints(data,constraints):
            nztot, row_block_heights = evaluate_blocking(grouping, nzcount, data["col_block_size"])
            data["padding"] = nztot - data["nonzeros"]
            if data["blocking_algo"] == 1: #structured
                data["padding"] /= 2
            data["avg_height"] = np.average(row_block_heights)
            datapoints.append(data)
            #print(data)
    return datapoints


def check_unique(l):
    return len(l) == len(set(l))

def add_curve(folder, x_name = "tau", y_name = "nonzeros_padding", constraints = {}, label = ""):
    datapoints = get_data_line(folder, constraints)
    x = [data[x_name] for data in datapoints]
    y = [data[y_name] for data in datapoints]
    L = sorted(zip(x,y), key=operator.itemgetter(0))
    x_s, y_s = zip(*L)
    if (not check_unique(x_s)):
        print("WARNING: plotting elements are not unique. Define better constraints")
    plt.plot(x_s, y_s, label = label, marker = "o")


algos = {}
algos["basic"] = {"blocking_algo": 0, "reorder": 0, "use_pattern": 0}
algos["pattern"] = {"blocking_algo": 0, "reorder": 0, "use_pattern": 1}
algos["structured"] = {"blocking_algo": 1, "reorder": 0, "use_pattern": 1}
algos["fixed"] = {"blocking_algo": 2, "reorder": 0}




for matrix_name in ["social","ia","soc-pocket","twitter"]:
    try:
        for block_size in [16,64,256]:
            x_name = "padding"
            y_name = "avg_height"

            print(f"Making {x_name} vs {y_name} image for graph {matrix_name}; block size = {block_size}")

            plt.figure()
            for algo, parameters in algos.items():
                parameters["col_block_size"] = block_size
                add_curve(f"results/{matrix_name}", x_name = x_name, y_name = y_name, constraints = parameters, label = algo)
            savename = f"images/reordering_curves_mat_{matrix_name}_X_{x_name}_Y_{y_name}_b_{block_size}.pdf"
            plt.xlabel(x_name)
            plt.ylabel(y_name)
            plt.yscale("log")
            plt.xscale("log")
            plt.legend()
            plt.savefig(savename,  bbox_inches='tight')
    except:
        print(f"could not create image for {matrix_name}")

