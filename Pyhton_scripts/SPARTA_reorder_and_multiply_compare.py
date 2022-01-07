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
    csvs["SA"] = "../results/test_cublas_reordering-saad-01-07-2022.csv";
    csvs["US"] = "../results/test_cublas_reordering-12-16-2021.csv";
    csvs["non_hierarchic"] = "../results/test_cublas_reordering-non-hierarchic-01-07-2022.csv";
    csvs["scalar"] = "../results/test_cublas_reordering-scalar-12-13-2021.csv";


    output_dir = "../images/";
                    
    
    dfs = {}
    for name in ["SA","US"]:
        dfs[name] = import_results(csvs[name])
        
    
    v_dict = {
            "US":     {"input_block_size" : 64, "merge_limit" : 0},
            "SA":     {"input_block_size" : 64, "merge_limit" : 0}   
            }
     
    compare_blocking_curve(dfs, v_dict, name =  "blocking_curve_compare_input_entries_US_SA")


    dfs = {}
    for name in ["scalar","US"]:
        dfs[name] = import_results(csvs[name])
        
    v_dict = {
            "US":     {"input_block_size" : 64, "merge_limit" : 0},
            "scalar":     {"input_block_size" : 64, "merge_limit" : 0}   
            }
     
    compare_blocking_curve(dfs, v_dict, name =  "blocking_curve_compare_input_entries_US_scalar")
    
    
    dfs = {}
    for name in ["SA","US"]:
        dfs[name] = import_results(csvs[name])
    vardict = {"input_block_size" : 64, "merge_limit" : 0}

    compare_heatmap(dfs["US"], dfs["SA"], vardict, vardict)

    
    