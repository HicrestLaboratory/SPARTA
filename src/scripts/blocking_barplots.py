# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import operator
import pandas as pd

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
    

def extract_data(filename):
    #header:matrix,rows,cols,nonzeros,blocking_algo,tau,row_block_size,col_block_size,use_pattern,sim_use_groups,sim_measure,reorder,exp_name,b_cols,warmup,exp_repetitions,multiplication_algo,n_streams,time_to_block,time_to_merge,time_to_compare,VBR_nzcount,VBR_nzblocks_count,VBR_average_height,merge_counter,comparison_counter,average_merge_tau,average_row_distance,avg_time_multiply,std_time_multiply,

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
        return data
    

def check_constraints(data_dict, constraints_dict):
    for constraint in constraints_dict:
        if not constraint in data_dict:
            print(f"WARNING: checking non existent property: {constraint} ")
            return False
        else:
            if constraints_dict[constraint] != data_dict[constraint]:
                return False 
    return True


def get_results(folder, constraints, variable):
    graph_names = []
    values = []
    for matrix_folder in glob.glob(f"{folder}/*"):
        graph_name = matrix_folder.split("/")[-1]
        graph_names.append(graph_name)
        valid_data = []
        for experiment in glob.glob(f"{matrix_folder}/*.txt"):
            data = extract_data(experiment)
            if check_constraints(data,constraints):
                valid_data.append(data[variable])
        values.append(min(valid_data))
    return graph_names, values



folder = "results/minitest"
savename = f"{folder}/../minitest_barplot.png"


variable = "VBR_nzcount"
plt.figure()
plt.xlabel("graphs")
plt.ylabel("# of nonzero 16 x 16 blocks")


barsize = 0.8
barpos = -barsize/2
increment = barsize/4
width = increment*0.9


hatches = {2 : "///", 5 : ".."}

for block_size in [16,32]:
    for algo, algoname in zip([2,5],("no-reordering","our reordering")):
        constraints = {"row_block_size" : block_size, "col_block_size": block_size, "blocking_algo" : algo}
        graph_names, values = get_results(folder, constraints, variable)
        print(block_size, algo, graph_names, values)
        x_pos = np.arange(barpos,len(graph_names) + barpos)
        print(x_pos)
        barpos += increment
        plt.bar(x_pos,values,label=f"{algoname}, block size = {block_size}", width = width, hatch = hatches[algo])
plt.legend()
plt.xticks(range(len(graph_names)), graph_names,rotation=90)
plt.savefig(savename,  bbox_inches='tight', dpi = 300)