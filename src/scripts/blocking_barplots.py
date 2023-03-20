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
import os as os

def get_dataframe(folder):
    df = pd.DataFrame()
    for matrix_folder in glob.glob(f"{folder}/*"):
        for experiment in glob.glob(f"{matrix_folder}/*"):
            try:
                df1 = pd.read_csv(experiment, skiprows = 0)
                df = pd.concat([df, df1], ignore_index=True, sort=False)
            except:
                print(f"CANNOT LOAD {experiment}")
    if df.empty:
        print("ALTERNATIVE LOAD", folder)
        for experiment in glob.glob(f"{folder}/*"):
            print(experiment)
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

def get_best_blockings(df, variable, variable_2, constraints = {}):
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


folder_name = "suitsparse_collection_2"
folder = f"results/{folder_name}"


variable = "block_density"
variable_2 = "tau"
ylabel = "Density"
savename = f"results/images/suitsparse_collection/suitsparse_collection_{variable}"
df = pd.DataFrame()
for fold_num in (3,4):
    folder_name = "suitsparse_collection_" + str(fold_num)
    folder = f"results/{folder_name}"
    df = pd.concat([df,get_dataframe(folder)],ignore_index=True, sort=False)

constraints = {}
for row_block_size in sorted(df["row_block_size"].unique()):
    for col_block_size in sorted(df["col_block_size"].unique()):
        plt.figure()
        plt.xlabel("graphs")
        plt.ylabel(ylabel)
        
        var_values = {};
        bars = 3
        barpos = -0.45
        increment = 0.9/bars
        width = increment*0.95

        for algo,algoname in zip((2,5),("blocking, no-reordering","blocking, reordering")):
            constraints["col_block_size"] = col_block_size
            constraints["row_block_size"] = row_block_size
            constraints["blocking_algo"] = algo
            res_df = get_best_blockings(df, variable = "VBR_nzblocks_count", variable_2 = "tau", constraints = constraints)
            matrices = [val.split("/")[-1].split(".")[0] for val in res_df["matrix"].values]
            taus = res_df["tau"].values
            res_df["density"] = res_df["nonzeros"].values/(res_df["rows"].values * res_df["cols"].values)
            res_df["block_density"] = res_df["nonzeros"].values/res_df["VBR_nzcount"].values
            var_values[algo] = res_df[variable].values
            print(algo, row_block_size, col_block_size, matrices, taus, var_values)
            x_pos = np.arange(barpos,len(matrices) + barpos)
            barpos += increment
            plt.bar(x_pos,var_values[algo],label=f"{algoname} ", width = width, hatch = "//", edgecolor = "black")
        if len(var_values) == 0: 
            plt.close()    
            continue
        
        dense_amp = np.nanmean(var_values[5]/var_values[2])
        dense_amp_orig = np.nanmean(var_values[5]/res_df["density"].values)
        print("AVG AMPLIFICATION", dense_amp)
        print("AVG AMPLIFICATION, original", dense_amp_orig)
        plt.text(1, 1, f"dense-amp: {dense_amp}, {dense_amp_orig}", bbox=dict(fill=False, edgecolor='red', linewidth=2))

        x_pos = np.arange(barpos,len(matrices) + barpos)
        barpos += increment
        plt.bar(x_pos, res_df["density"].values, label = "no blocking (original)", width = width, hatch = "..", color = "white", edgecolor='black')
        

        plt.legend()
        plt.yscale("log")
        plt.title(f"block height = {row_block_size}, block width = {col_block_size}")
        plt.xticks(range(len(matrices)), matrices, rotation=90)
        plt.savefig(savename + f"_{row_block_size}_{col_block_size}.png",  bbox_inches='tight', dpi = 300)
        plt.close()    