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
    

def extract_data(filename, skip = 0):
    #header:matrix,rows,cols,nonzeros,blocking_algo,tau,row_block_size,col_block_size,use_pattern,sim_use_groups,sim_measure,reorder,exp_name,b_cols,warmup,exp_repetitions,multiplication_algo,n_streams,time_to_block,time_to_merge,time_to_compare,VBR_nzcount,VBR_nzblocks_count,VBR_average_height,merge_counter,comparison_counter,average_merge_tau,average_row_distance,avg_time_multiply,std_time_multiply,

    with open(filename) as infile:
        complete = False
        for idx, line in enumerate(infile):
            line = line.split(",")[:-1]
            if idx < skip:
                continue
            if idx == skip:
                vars = line
            elif idx == skip+1:
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
    
    best_values = {}
    relevant_variable = "VBR_longest_row"

    for matrix_folder in glob.glob(f"{folder}/*"):
        graph_name = matrix_folder.split("/")[-1]
        valid_data = []

        best_values[graph_name] = []
        current_min = float('inf')
        try:
            for experiment in glob.glob(f"{matrix_folder}/*.txt"):
                data = extract_data(experiment)
                data["effective_density"] = data["nonzeros"]/data["VBR_nzcount"]
                data["true_density"] = data["nonzeros"]/(data["rows"]*data["cols"])
                data["relative_density"] = data["effective_density"]/(data["true_density"])
                data["nz_per_block"] = data["nonzeros"]/data["VBR_nzblocks_count"]

                if check_constraints(data,constraints):
                    valid_data.append(data[variable])

                    if constraints["blocking_algo"] == 5:
                        if (data[relevant_variable] < current_min):
                            current_min = data[relevant_variable]
                            best_values[graph_name] = data["tau"]
                        
            values.append(max(valid_data))
            graph_names.append(graph_name)
        except: 
            print(f"ERROR LOADING DATA FOR GRAPH {graph_name}", constraints)

    print("FOUND BEST", relevant_variable, " : ", constraints["col_block_size"], best_values)
    return graph_names, np.array(values)




folder = "results/minitest"
savename = f"{folder}/../denseAMP_miniset_relative"



variable = "effective_density"
ylabel = "density amplification relative to natural blocking"





def get_dataframe(folder):
    df = pd.DataFrame()
    for matrix_folder in glob.glob(f"{folder}/*"):
        for experiment in glob.glob(f"{matrix_folder}/*"):
            try:
                df1 = pd.read_csv(experiment, skiprows = 0)
                df = pd.concat([df, df1], ignore_index=True, sort=False)
            except:
                print(f"CANNOT LOAD {experiment}")
    return df        

def get_line(df, constraints):
    query = ""
    for key in constraints:
        query += f"{key} == {constraints[key]} &"
    query = query[:-2]
    return df.query(query)

def apply_function_per_matrix(df, variable, variable_2, constraints = {}):
    tmp_df = get_line(df, constraints)
    result_df = pd.DataFrame()
    for matrix_path in df["matrix"].unique():
        matrix_name = matrix_path.split("/")[-1]
        mat_df = tmp_df[tmp_df["matrix"] == matrix_path]
        mat_df = mat_df[mat_df[variable] == mat_df[variable].min()]
        mat_df = (mat_df[mat_df[variable_2] == mat_df[variable_2].min()])
        try:
            result_df = pd.concat([result_df, mat_df], ignore_index=True, sort=False)
        except:
            print("found no values for matrix ", matrix_name)
    return(result_df)


folder = "results/minitest_hamming"
variable = "VBR_longest_row"
variable_2 = "tau"
savename = f"results/minitest_hamming_{variable}"
df = get_dataframe(folder)
constraints = {}
for block_size in [16,32,64,128]:
    plt.figure()
    plt.xlabel("graphs")
    plt.ylabel(ylabel)
    barsize = 0.4
    barpos = -barsize/2
    increment = barsize
    width = increment*0.9
    for algo,algoname in zip((2,5),("no-reordering","reordering")):
        constraints["col_block_size"] = block_size
        constraints["blocking_algo"] = algo
        res_df = apply_function_per_matrix(df, variable = variable, variable_2 = variable_2, constraints = constraints)
        matrices = res_df["matrix"].values
        taus = res_df["tau"].values
        var_values = res_df[variable].values
        print(algo, block_size, matrices, taus, var_values)
        x_pos = np.arange(barpos,len(matrices) + barpos)
        barpos += increment
        plt.bar(x_pos,var_values,label=f"{algoname}, block size = {block_size}", width = width)
    plt.legend()
    plt.xticks(range(len(matrices)), matrices, rotation=90)
    plt.savefig(savename + f"_{block_size}.png",  bbox_inches='tight', dpi = 300)