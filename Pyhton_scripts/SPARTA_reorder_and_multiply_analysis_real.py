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
                    
results_df = import_results(input_csv)

results_df["input_source"] = results_df.apply(lambda x: x['input_source'].split("/")[-1], axis = 1)

plt.rcParams['font.size'] = 14

do_all = False
if do_all:
    
    ignore = ["input_source","rows","cols"];
    fixed = {"B_cols" : 4096, "similarity_func" : "'jaccard'", "reorder_algorithm": "'saad_blocks'", "hierarchic_merge" : 1, "merge_limit" : 0};
    for values in generate_exp_iterator(results_df, ignore = ignore, fixed = fixed):
        variables_dict = dict(zip(experimental_variables, list(values)))
        print(variables_dict)
        #try:
        bar_plot(results_df, variables_dict);
        #except:
    
    
    ignore = ["cols","B_cols","epsilon"];
    fixed = {"similarity_func" : "'jaccard'", "reorder_algorithm": "'saad_blocks'"};
    for values in generate_exp_iterator(ignore = ignore, fixed = fixed):
        variables_dict = dict(zip(experimental_variables, list(values)))
        #try:
        reorder_curve(variables_dict);
        #except:

save_folder = "../images/paper_images_real/"

variables_dict = {"reorder_algorithm": "'saad_blocks'", "merge_limit" : -1, "epsilon" : 0.5, "B_cols" : 4096, "hierarchic_merge": 1}
#real_blocking_curve(results_df, variables_dict, variable = "input_source")        
#bar_plot_together(results_df, variables_dict, save_folder = save_folder, name = "very_small");

for graph in results_df["input_source"].unique():
    variables_dict = {"input_source": "'" + graph + "'" , "reorder_algorithm": "'saad_blocks'", "merge_limit" : -1, "epsilon" : 0.5, "B_cols" : 4096, "hierarchic_merge": 1}
    q = build_query(variables_dict);
    this_df = results_df.query(q)
    print(graph.split(".")[0], "&")
    print(this_df["rows"].values[0], "&")
    print(this_df["total_nonzeros"].values[0], "&")
    print(format(this_df["input_density"].values[0],".4f"), "\\\\")



