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

                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_cublas_reordering-12-16-2021.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;
    results_df = import_results(input_csv)

    do_all_images = False;
    
    plt.rcParams['font.size'] = 15

    
    if do_all_images:
        ignore = ["input_entries_density","input_blocks_density", "epsilon", "input_block_size"];
        fixed = {"similarity_func" : "'jaccard'", "reorder_algorithm": "'saad_blocks'"};
        for values in generate_exp_iterator(ignore = ignore, fixed = fixed):
            variables_dict = dict(zip(experimental_variables, list(values)))
            try:
                performance_heatmap(results_df, variables_dict);
                epsilon_heatmap(results_df, variables_dict)
                reorder_heatmap(results_df, variables_dict)
                delta_heatmap(results_df, variables_dict)
            except Exception as e:
                print(e, variables_dict)
            
            
    #paper images
    
    
    
    save_folder = "../images/paper_images/"
    for B_cols in [2048, 4096, 8192]:
        variables_dict = {"rows": 8192, "cols": 8192, "B_cols": B_cols, "input_block_size" : 64, "algo_block_size" : 64, "reorder_algorithm": "'saad_blocks'", "merge_limit" : -1}
        performance_heatmap(results_df, variables_dict, save_folder = save_folder);
        epsilon_heatmap(results_df, variables_dict, save_folder = save_folder)
    
    variables_dict = {"rows": 8192, "cols": 8192, "B_cols": 8192, "input_block_size" : 64, "algo_block_size" : 64, "reorder_algorithm": "'saad_blocks'", "merge_limit" : 0}
    reorder_heatmap(results_df, variables_dict, save_folder = save_folder)
    delta_heatmap(results_df, variables_dict, save_folder = save_folder)

    
    variables_dict = {"rows": 8192, "cols": 8192, "B_cols": 8192, "input_block_size" : 64, "algo_block_size" : 64,"input_blocks_density" : 0.1, "reorder_algorithm": "'saad_blocks'", "merge_limit" : 0}
    blocking_curve(results_df, variables_dict, save_folder = save_folder)

    variables_dict = {"rows": 8192, "cols": 8192, "B_cols": 8192, "input_block_size" : 64, "algo_block_size" : 64,"input_blocks_density" : 0.1, "input_entries_density" : 0.1, "reorder_algorithm": "'saad_blocks'"}
    blocking_curve(results_df, variables_dict, variable = "merge_limit", save_folder = save_folder, labels = ["Theory", "None", "2.", "4."], title = "original in-block density = 0.1", savename = "merge_comparison")