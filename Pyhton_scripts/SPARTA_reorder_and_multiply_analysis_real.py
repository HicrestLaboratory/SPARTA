# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:16:49 2021

@author: Paolo
"""


import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns
import itertools as itr
from scipy import interpolate
from analysis_tools import *

                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_cublas_reordering-real-verysmall-01-11-2022.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;


"""                    
dfs = {}
input_csv = "../results/test_cublas_reordering-real-verysmall-01-11-2022.csv"
this_df = import_results(input_csv)
this_df = this_df[this_df["rows"] < 15000]
dfs["very_small"] = this_df
input_csv = "../results/test_cublas_reordering-real-small-01-13-2022.csv"
this_df = import_results(input_csv)
this_df = this_df[this_df["rows"] != 21608]
dfs["small"] = this_df[this_df["rows"] > 4000]
"""

dfs = {}
input_csv = "../results/test_cublas_reordering-real-verysmall-01-11-2022.csv"
this_df_1 = import_results(input_csv)
input_csv = "../results/test_cublas_reordering-real-small-01-13-2022.csv"
this_df_2 = import_results(input_csv)
this_df_merge = pd.concat([this_df_1,this_df_2])
this_df = this_df_merge.sort_values("rows", ascending = True)



dfs["small"] = this_df[this_df["rows"] < 4000]
dfs["medium"] = this_df[this_df["rows"] > 4000]
dfs["medium"] = dfs["medium"][dfs["medium"]["rows"] < 20000]
dfs["big"] = this_df[this_df["rows"] > 20000]



#dfs["small_plus"] = this_df[this_df["rows"] > 15000]

#dfs["small"].replace("movielens-10m-noRatings.el","movielens-10-nR.el",inplace = True)

for df_name,df in dfs.items():
    df["input_source"] = df.apply(lambda x: x['input_source'].split("/")[-1], axis = 1)

plt.rcParams['font.size'] = 14

save_folder = "../images/paper_images_real/"

variables_dict = {"reorder_algorithm": "'saad_blocks'", "merge_limit" : 0, "epsilon" : 0.1, "B_cols" : 4096, "hierarchic_merge": 1}
#real_blocking_curve(dfs["very_smalle"], variables_dict, variable = "input_source")        
#bar_plot_together(dfs["very_small"], variables_dict, save_folder = save_folder, name = "very_small");


for df_name,df in dfs.items():
    bar_plot_together(df, variables_dict, save_folder = save_folder, name = df_name);

for df_name,df in dfs.items():
    df.sort_values("rows", ascending = True, inplace = True)

    for graph in df["input_source"].unique():
        variables_dict = {"input_source": "'" + graph + "'" , "reorder_algorithm": "'saad_blocks'", "merge_limit" : 0, "epsilon" : 0.1, "B_cols" : 4096, "hierarchic_merge": 1}
        q = build_query(variables_dict);
        this_df = df.query(q)
        print(graph.split(".")[0], "&")
        
        print(this_df["rows"].values[0], "&")
        print(this_df["total_nonzeros"].values[0], "&")
        print(format(this_df["input_density"].values[0]*100,".3f") + "\\%", "\\\\")



