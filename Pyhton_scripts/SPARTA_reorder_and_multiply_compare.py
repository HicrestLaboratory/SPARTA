# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 17:28:51 2022

@author: Paolo
"""


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
    
    parser.add_argument("--input-csv", default="../results/test_cublas_reordering-saad-01-07-2022.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--original-csv", default="../results/test_cublas_reordering-12-16-2021.csv",
        help="file that contains the already run experiments for the original")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    csvs = {};
    csvs["SA"] = "../results/test_cublas_reordering-saad-01-10-2022.csv";
    csvs["1-SA"] = "../results/test_cublas_reordering-synthetic-01-11-2022.csv";
    csvs["non_hierarchic"] = "../results/test_cublas_reordering-non-hierarchic-01-11-2022.csv";
    csvs["cosine"] = "../results/test_cublas_reordering-scalar-01-10-2022.csv";
    
    
    
    save_folder = "../images/paper_images_new/"
    output_dir = "../images/";
                    
    plt.rcParams['font.size'] = 16

    dfs = {}
    for name in csvs:
        dfs[name] = import_results(csvs[name])
    
    
    v_dict = {}
    v_dict["1-SA"] = {"input_block_size" : 64, "merge_limit" : 0, "B_cols" : 2048}           
    v_dict["SA"] = {"input_block_size" : 64, "merge_limit" : 0}   

    
    
    
    compare_blocking_points(dfs["SA"], "SA", v_dict["SA"], dfs["1-SA"], "1-SA", v_dict["1-SA"], name =  "blocking_points_compare_input_entries_US_SA", save_folder = save_folder)
    
    
    
    v_dict["cosine"] = {"input_block_size" : 64, "merge_limit" : 0}   

    compare_blocking_points(dfs["1-SA"], "jaccard", v_dict["1-SA"], dfs["cosine"], "cosine", v_dict["cosine"], name =  "blocking_curve_compare_input_entries_US_cosine", save_folder = save_folder)


    dfs = {}
    for name in ["SA","1-SA"]:
        dfs[name] = import_results(csvs[name])
        
    v_dict = {
            
            "SA":     {"input_block_size" : 64, "merge_limit" : 0, "input_blocks_density": 0.1, "B_cols": 128},
            "1-SA":     {"input_block_size" : 64, "merge_limit" : 0, "B_cols" : 2048}
            }
     
    
    
    compare_blocking_curves(dfs, 0.1 ,[0.01,0.1, 0.2, 0.5], v_dict, name =  "blocking_curve_compare_US_SA", save_folder = save_folder)
    
    blocking_curve(dfs["SA"], v_dict["SA"], values = [0.01,0.02,0.1,0.2,0.5], name = "saad_curve_motivation",  save_folder = save_folder)
    
    
    dfs = {}
    for name in csvs:
        dfs[name] = import_results(csvs[name])
    vardict_us = {"input_block_size" : 64, "merge_limit" : -1}
    vardict_them = {"input_block_size" : 64, "merge_limit" : 0}

    compare_heatmap(dfs["1-SA"], dfs["cosine"], vardict_us, vardict_them, save_folder = save_folder)
    compare_heatmap(dfs["1-SA"], dfs["SA"], vardict_us, vardict_them, save_folder = save_folder)

    
    variables_dict = {"rows": 8192, "cols": 8192, "input_blocks_density" : 0.1, "merge_limit" : 0}
    blocking_curve(dfs["SA"], variables_dict, save_folder = "../images/tests")
    
    variables_dict = {"rows": 8192, "cols": 8192, "B_cols" : 8192, "input_blocks_density" : 0.1, "merge_limit" : 0}
    blocking_curve(dfs["1-SA"], variables_dict, save_folder = "../images/tests")