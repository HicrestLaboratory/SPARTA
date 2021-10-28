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
                 
   
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_cublas_results_ultrasparse_19_10.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;
                    
    
    
input_csv = "../results/a100_test01.csv"
results_df = pd.read_csv(input_csv);

columns = ['exp_name',
'input_type',
'input_source',
'rows',
'cols',
'B_cols',
'B_density',
'total_nonzeros',
'input_blocks_density',
'input_entries_density',
'input_block_size',
'reorder_algorithm',
'algo_block_size',
'epsilon',
'similarity_func',
'hierarchic_merge',
'scramble',
'Warmup',
'Repetitions',
'Algorithm',
'VBSmm_mean(ms)',
'VBSmm_std',
'cusparse_spmm_mean(ms)',
'cusparse_spmm_std']

numerics = ["rows",'cols',
'B_cols',
'B_density',
'total_nonzeros',
'input_blocks_density',
'input_entries_density',
'input_block_size',
'algo_block_size',
'epsilon',
'VBSmm_mean(ms)',
'VBSmm_std',
'cusparse_spmm_mean(ms)',
'cusparse_spmm_std']

results_df[numerics] = results_df[numerics].apply(pd.to_numeric)

results_df["density"] = results_df.apply(lambda x: x['input_blocks_density']*x["input_entries_density"], axis=1)
results_df["sp_vs_cu"] = results_df.apply(lambda x: x['cusparse_spmm_mean(ms)']/x["VBSmm_mean(ms)"], axis=1)
results_df["combined_std"] = results_df.apply(lambda x: x["sp_vs_cu"]/((x["VBSmm_mean(ms)"]/(3*x['VBSmm_std']))**2 + (x["cusparse_spmm_mean(ms)"]/(3*x['cusparse_spmm_std']))**2)**(0.5), axis=1)
        
results_df["in-block-density"] = results_df.apply(lambda x: x['density']/x["input_blocks_density"], axis=1)
def winner(x):
    win = 1;
    if (x["VBSmm_mean(ms)"] < x["cusparse_spmm_mean(ms)"]):
        win = 0;
    else:
        win = 2;
    return win;

results_df["winner"] = results_df.apply(winner, axis = 1);

print(results_df["input_entries_density"].unique())
print(results_df["input_blocks_density"].unique())


to_display = ["rows","cols", 
              "input_entries_density",
              "input_blocks_density",
              "rows",
              "cols",
              "B_cols",
              "input_block_size",
              ]
for var in to_display:
    print(var, results_df[var].unique());

def make_heatmap(cols, rows, b_cols, block_size, save_folder):
    q = "cols == {} and rows == {} and B_cols == {} and input_block_size == {}".format(cols, rows, b_cols, block_size)
    heatmap_df = results_df.query(q)[["input_entries_density","input_blocks_density","sp_vs_cu"]]
    heatmap_df = heatmap_df.set_index(['input_entries_density', 'input_blocks_density']).sp_vs_cu.unstack(0)
    heatmap_df.sort_index(level=0, ascending=False, inplace=True)
    
    sns.set(font_scale=1)
    cmap = sns.diverging_palette(0,255,sep=1, as_cmap=True)
    
    ax = sns.heatmap(heatmap_df, annot=True, linewidths=.5, cbar_kws={'label': 'speed-up vs cuSparse'}, cmap = cmap, center = 1, vmin = 0, vmax = 5)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    
    plt.xlabel("Density inside nonzero blocks") 
    plt.ylabel("Density of nonzero blocks") 
    
    plt.title("M,K,N = {},{},{} \n block_size = {}".format(cols,rows,b_cols,block_size))
    
    savename = save_folder + "landscape_heatmap_r{}_c{}_k{}_b{}.jpg".format(rows, cols, b_cols, block_size);
    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    plt.close()

for cols in [2048,4096,8192]:
    for rows in [2048,4096,8192]:
        for b_cols in [2048,4096,8192,16384]:
            for block_size in [32,64,128]:
                try:
                    make_heatmap(cols,rows,b_cols,block_size, save_folder = "../images/performance_landscape/");
                except:
                    print("could not make image for cols = {}, rows = {}, B_cols = {}, block_size = {}".format(cols,rows,b_cols,block_size))






def speedup_curve(rows,cols,B_cols,b_size,b_density, save_folder = "../images/reoder_curve/"):
        
    plt.figure()
    ax = plt.gca();
    
    y = "sp_vs_cu";
    x = "input_entries_density"
    q = "rows == {} and cols == {} and B_cols == {} and input_block_size == {} and input_blocks_density == {}".format(rows, cols, B_cols, b_size, b_density)
    results_df.query(q).plot(x = x, y = y, color = "red", ax = ax);
    ax.set(ylabel='speedup vs cusparse', xlabel='density inside blocks',
           title='Comparison with cusparse \n N, K, M = {}, {}, {} \n block size = {} \n input_blocks_density = {}'.format(rows, cols, B_cols, b_size, b_density))
    ax.axhline(1)
    x = results_df.query(q)["input_entries_density"].tolist();
    y = np.array(results_df.query(q)["sp_vs_cu"].tolist());
    errors = np.array(results_df.query(q)["VBSmm_std"].tolist());
    upper = y + errors;
    lower = y - errors;

    ax.set_xscale("log")
    ax.set_xlim([0.0008, 0.31]);
    ax.set_ylim([0, 10]);

    
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
    
                    
    ax.fill_between(x, lower, upper, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    
    plt.title("M,K,N = {},{},{} \n block size = {} \n fraction of nonzero blocks = {}".format(rows,cols,B_cols,b_size,b_density))
    
    savename = save_folder + "reorder_curve_r{}_c{}_k{}_bs{}_bd{}.jpg".format(rows,cols,B_cols,b_size,b_density);                
    
    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    
    plt.show()
    plt.close()


for cols in [8192,]:
    for rows in [8192,]:
        for b_cols in [2048, 4096, 8192, 16384]:
            for block_size in [32,64,128]:
                for b_density in [0.05, 0.2, 0.4]:
                    try:
                        speedup_curve(rows,cols,b_cols,block_size,b_density, save_folder = "../images/reorder_curve/");
                    except:
                        print("could not make image for cols = {}, rows = {}, B_cols = {}, block_size = {}".format(cols,rows,b_cols,block_size))
