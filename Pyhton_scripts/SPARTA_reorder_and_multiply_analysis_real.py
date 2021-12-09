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
    
    parser.add_argument("--input-csv", default="../results/test_cublas_reordering-real-12-07-2021.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;
                    
results_df = pd.read_csv(input_csv);

columns = {'exp_name' : "experiment name",
           'input_type': "type of input", 
           'input_source': "input matrix", 
           'rows' : "K", 
           'cols': "M", 
           'B_cols': "N",
       'B_density': "density of B", 
       'total_nonzeros': "total nonzeros", 
       'input_blocks_density': "fraction of nonzero blocks",
       'input_entries_density': "density inside nonzero blocks", 
       'input_block_size': "block size", 
       'reorder_algorithm': "reorder algorithm",
       'algo_block_size': "block size", 
       'epsilon' : "tau", 
       'similarity_func': "similarity function", 
       'hierarchic_merge': "hierarchic merge",
       'merge_limit': "merge limit", 
       'scramble': "scramble", 
       'Warmup': "warmup experiments", 
       'Repetitions': "number of repetitions", 
       'Algorithm': "algo sel",
       'VBS_avg_nzblock_height': "avg. block height in the reordered matrix",
       'VBS_avg_nzblock_height_error': "error of the avg block height",
       'VBS_total_nonzeros': "total nonzero area in the reordered matrix",
       'VBS_total_nonzeros_error': "total nonzero area error", 
       'VBS_block_rows': "block rows",
       'VBS_block_rows_error': "block rows error", 
       'VBS_nz_blocks': "number of nonzero blocks", 
       'VBS_nz_blocks_error': "error on the number of nonzero blocks",
       'avg_skipped': "avg number of skipped comparisons", 
       'skipped_std': "error of skipped", 
       'avg_comparisons': "avg number of row comparisons", 
       'comparisons_std': "error in row comparisons",
       'VBSmm_algo_mean(ms)': "time to complete a block multiplication", 
       'VBSmm_algo_std': "error in block multiplication", 
       'cusparse_spmm_mean(ms)': "time to complete the cusparse multiplication",
       'cusparse_spmm_std': "error in cusparse multiplication"
       }

numerics = [
           'rows', 
           'cols', 
           'B_cols',
       'B_density',
       'total_nonzeros', 
       'input_blocks_density',
       'input_entries_density', 
       'input_block_size', 
       'algo_block_size', 
       'epsilon', 
       'Warmup', 
       'Repetitions',
       'merge_limit',
       'Algorithm',
       'VBS_avg_nzblock_height', 
       'VBS_avg_nzblock_height_error',
       'VBS_total_nonzeros', 
       'VBS_total_nonzeros_error', 
       'VBS_block_rows',
       'VBS_block_rows_error', 
       'avg_skipped', 
       'skipped_std', 
       'avg_comparisons', 
       'comparisons_std',
       'VBS_nz_blocks', 
       'VBS_nz_blocks_error',
       'VBSmm_algo_mean(ms)',
       'VBSmm_algo_std',
       'cusparse_spmm_mean(ms)',
       'cusparse_spmm_std']

results_df[numerics] = results_df[numerics].apply(pd.to_numeric)

results_df["input_density"] = results_df.apply(lambda x: x['total_nonzeros']/(x["rows"]*x["cols"]), axis=1)

results_df["output_in_block_density"] = results_df.apply(lambda x: x['total_nonzeros']/x["VBS_total_nonzeros"], axis=1)

results_df["relative_density"] = results_df.apply(lambda x: x['output_in_block_density']/x["input_entries_density"], axis=1)

results_df["true_relative_density"] = results_df.apply(lambda x: x['output_in_block_density']/x["input_density"], axis=1)

results_df["sp_vs_cu"] = results_df.apply(lambda x: x['cusparse_spmm_mean(ms)']/x["VBSmm_algo_mean(ms)"], axis=1)

#results_df["input_source"] = results_df.apply(lambda x: x['input_source'].split("/")[-1], axis = 1)


experimental_variables = [
              "input_entries_density",
              "input_blocks_density",
              "rows",
              "cols",
              "B_cols",
              "input_block_size",
              "algo_block_size",
              "epsilon",
              "similarity_func",
              "scramble",
              "hierarchic_merge",
              "merge_limit",
              "reorder_algorithm"
        ]

for var in experimental_variables:
    print(var, results_df[var].unique());

def build_query(fixed):
    #build a query from the fixed dictionary
    q = ""
    for k, val in fixed.items():
        if val != "any":
            q += k + " == " + str(val) + " and ";
    q = q[0:-5];
    return q;   




fixed = {}
for var in experimental_variables:
    fixed[var] = "any";


def generate_exp_iterator(ignore = [], fixed = {}):
    value_lists = []
    for var in experimental_variables:
        if (var in fixed.keys()):
            value_lists.append([fixed[var],])
        elif (var in ignore):
            value_lists.append(["any",])
        else:
            value_lists.append(results_df[var].unique())
            
    return itr.product(*value_lists);


def make_title(variables_dict, ignore = ["algo_block_size","scramble","reorder_algorithm"]):
    q = ""
    for k, val in variables_dict.items():
        if val != "any" and k not in ignore: 
            q += columns[k] + " = " + str(val) + " \n ";
    return q;

def add_to_query(var, val):
    return " and " + var + "==" + str(val);

def make_savename(name, variables_dict, ignore = ["input_block_size",]):
    q = name
    for k, val in variables_dict.items():
        if val != "any" and k not in ignore: 
            q += "_" + k[0:2] + str(val);
    return q + ".jpg";


def bar_plot(variables_dict, save_folder = "../images/performance_real/", name = "real_barplots"):
    
    if variables_dict["reorder_algorithm"] == "saad": 
            name += "_saad";
    
    graphs = results_df["input_source"].unique();
    measures = ['cusparse_spmm_mean(ms)','VBSmm_algo_mean(ms)']
    
    plt.style.use('grayscale')    
    fig, ax = plt.subplots(1, figsize = (len(graphs), 4), sharex=True);


    colormap = plt.cm.Greys
    ax.tick_params(axis='x', colors='grey')
    ax.grid(b=True, which='major', color='#999999', linestyle='-', alpha=0.3, linewidth=0.6)
    ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.4)
    ax.set_prop_cycle(color=[colormap(i) for i in np.linspace(0, 1,len(measures))])

    space = 0.2;
    width = (1. - space)/len(measures);

    q = build_query(variables_dict)
    results_df.query(q).plot.bar(x = "input_source",  y = measures, width = width, edgecolor = "black", ax = ax);
     
    plt.ylabel("multiplication time (ms)");
    plt.xlabel("input graph");


    plt.title(make_title(variables_dict, ignore = ["scramble","input_source","input_block_size","input_blocks_density","input_entries_density"]))
    
    savename = make_savename(name,variables_dict)
    plt.savefig(save_folder + savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()


def reorder_curve(variables_dict, save_folder = "../images/performance_real/", name = "reorder_curve_real"):
   
    if variables_dict["reorder_algorithm"] == "saad": 
            name += "_saad";
            
    plt.style.use('grayscale')    

    q = build_query(variables_dict)
    results_df.query(q).plot(x = "VBS_avg_nzblock_height", y = "output_in_block_density", kind = "scatter");
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
    #plt.xlim(0,200)
    #plt.ylim(0,2)
    plt.axvline(variables_dict["algo_block_size"], alpha = 0.2, color = "red")
    plt.title(make_title(variables_dict, ignore = ["scramble","input_block_size","input_blocks_density","input_entries_density"]))
    
    plt.xlabel("Average height of nonzero blocks");
    plt.ylabel("Average density inside nonzero blocks \n (relative to original blocking)");
    
    savename = make_savename(name,variables_dict)

    #plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    plt.close()
    


ignore = ["input_source","rows","cols"];
fixed = {"similarity_func" : "'jaccard'", "reorder_algorithm": "'saad_blocks'"};
for values in generate_exp_iterator(ignore = ignore, fixed = fixed):
    variables_dict = dict(zip(experimental_variables, list(values)))
    #try:
    bar_plot(variables_dict);
    #except:


ignore = ["cols","B_cols","epsilon"];
fixed = {"similarity_func" : "'jaccard'", "reorder_algorithm": "'saad_blocks'"};
for values in generate_exp_iterator(ignore = ignore, fixed = fixed):
    variables_dict = dict(zip(experimental_variables, list(values)))
    #try:
    reorder_curve(variables_dict);
    #except: