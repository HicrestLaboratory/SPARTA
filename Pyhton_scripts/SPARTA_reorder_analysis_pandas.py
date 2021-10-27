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
from scipy import interpolate


saad = False
if saad:
    input_file = "../results/test_reordering_blocked_synth_saad_27_10.csv"
else:
    input_file = "../results/test_reordering_blocked_synth_27_10.csv"
#    input_file = "../results/real_reordering_results_25_10.csv"

                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_reordering_blocked_synth_2_saad.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;
                    
results_df = pd.read_csv(input_file);

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
       'VBS_avg_nzblock_height', 
       'VBS_avg_nzblock_height_error',
       'VBS_total_nonzeros', 
       'VBS_total_nonzeros_error', 
       'VBS_block_rows',
       'VBS_block_rows_error', 
       'VBS_nz_blocks', 
       'VBS_nz_blocks_error',
       'VBS_min_block_H', 
       'VBS_min_block_H_error', 
       'VBS_max_block_H',
       'VBS_max_block_H_error']

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
       'Algorithm',
       'VBS_avg_nzblock_height', 
       'VBS_avg_nzblock_height_error',
       'VBS_total_nonzeros', 
       'VBS_total_nonzeros_error', 
       'VBS_block_rows',
       'VBS_block_rows_error', 
       'VBS_nz_blocks', 
       'VBS_nz_blocks_error',
       'VBS_min_block_H', 
       'VBS_min_block_H_error', 
       'VBS_max_block_H',
       'VBS_max_block_H_error']

results_df[numerics] = results_df[numerics].apply(pd.to_numeric)

results_df["input_density"] = results_df.apply(lambda x: x['total_nonzeros']/(x["rows"]*x["cols"]), axis=1)

results_df["output_in_block_density"] = results_df.apply(lambda x: x['total_nonzeros']/x["VBS_total_nonzeros"], axis=1)

results_df["relative_density"] = results_df.apply(lambda x: x['output_in_block_density']/x["input_entries_density"], axis=1)

results_df["true_relative_density"] = results_df.apply(lambda x: x['output_in_block_density']/x["input_density"], axis=1)


results_df = results_df[results_df["input_blocks_density"] != 0.01]

to_display = ["input_entries_density",
              "input_blocks_density",
              "rows",
              "cols",
              "input_block_size",
              "algo_block_size",
              "epsilon",
              "similarity_func",
              "scramble",
              "hierarchic_merge"
        ]

for var in to_display:
    print(var, results_df[var].unique());

def build_query(fixed):
    #build a query from the fixed dictionary
    q = ""
    for k, val in fixed.items():
        if val != "any":
            q += k + " == " + str(val) + " and ";
    q = q[0:-5];
    return q;




def reorder_curve(rows = 2048, cols = 2048, block_size = 64, i_density = 0.1, b_density = 0.1, similarity = "'scalar'", save_folder = "../images/reorder_curve/", name = "reorder_curve", saad = False):
    fixed = {
        "input_entries_density": i_density,
        "input_blocks_density": b_density,
        "rows": str(rows),
        "cols": str(cols),
        "input_block_size": str(block_size),
        "algo_block_size": str(block_size),
        "epsilon": "any",
        "similarity_func": similarity,
        "scramble": "1"
    }
    
    if saad: 
        fixed["algo_block_size"] = 1;
        name = "reorder_curve_saad";

    fixed
    q = build_query(fixed)
    results_df.query(q).plot(x = "VBS_avg_nzblock_height", y = "relative_density", kind = "scatter");
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
    #plt.xlim(0,1000)
    #plt.ylim(0,2)
    plt.axvline(block_size, alpha = 0.2, color = "red")
    plt.axhline(1, alpha = 0.2, color = "red")
    plt.title("M,N = {},{}\n block_size = {} \n density inside nonzero blocks = {} \n fraction of nonzero-blocks = {} \n similarity function = {}".format(cols,rows,block_size,  i_density, b_density, similarity))
    plt.xlabel("Average height of nonzero blocks");
    plt.ylabel("Average density inside nonzero blocks \n (relative to original blocking)");
    
    savename = save_folder + name + "_r{}_c{}_b{}_d{}_bd{}_{}.jpg".format(cols,rows,block_size,  i_density, b_density, similarity);

    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    plt.close()
    
    
if saad:
    name = "saad_reorder_curve"    
else:
    name = "reorder_curve";
    
for cols in [2048,]:
    for rows in [2048,]:
        for block_size in [32,64]:
            for i_density in [0.01, 0.05, 0.1, 0.3]:
                for b_density in [0.1,0.2,0.3,0.4,0.5]:
                    for similarity in ["'scalar'","'jaccard'"]:
                        try:
                            reorder_curve(cols,rows,block_size, i_density, b_density, similarity, name = "reorder_curve", saad = saad);
                        except:
                            print("could not make image for cols = {}, rows = {}, block_size = {}".format(cols,rows,block_size))




def real_reorder_curve(graph = None, block_size = 64, similarity = "'scalar'", save_folder = "../images/reorder_curve/", name = "real_reorder_curve", saad = False):
    fixed = {
        "input_source": "'" + str(graph) + "'",
        "algo_block_size": str(block_size),
        "epsilon": "any",
        "similarity_func": similarity,
        "scramble": "1"
    }
    
    if saad: 
        fixed["algo_block_size"] = 1;
        name = "reorder_curve_saad";


    graph_name = graph.split("/")[-1].split(".")[0]
    fixed
    q = build_query(fixed)
    results_df.query(q).plot(x = "VBS_avg_nzblock_height", y = "true_relative_density", kind = "scatter");
    dens = str(results_df.query(q)["input_density"].unique()[0])[0:6];
    rows = str(results_df.query(q)["rows"].unique()[0]);
    cols = str(results_df.query(q)["cols"].unique()[0]);

    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
    plt.xlim(0,100)
    plt.ylim(0,20)
    plt.title("matrix: {}, \n M, K = {},{}, \n density = {} \n block size = {} \n".format(graph_name, rows, cols, dens, block_size, similarity))
    plt.xlabel("Average height of nonzero blocks");
    plt.ylabel("Average density inside nonzero blocks \n (relative to original density)");
    
    savename = save_folder + name + "{}_bs{}_{}.jpg".format(graph_name,block_size, similarity);

    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    plt.close()
    
    
for graph in results_df["input_source"].unique():
    for algo_block_size in [64,]:
        for sim in ["'scalar'",]:
            try:
                real_reorder_curve(graph,algo_block_size,sim);
            except: 
                print("could not process", graph)


def curve_comparison(rows = 2048, cols = 2048, block_size = 64, i_density = 0.1, b_density = 0.1, save_folder = "../images/reorder_curve/", name = "reorder_curve_comparison"):
    fixed = {
        "input_entries_density": i_density,
        "input_blocks_density": b_density,
        "rows": str(rows),
        "cols": str(cols),
        "input_block_size": str(block_size),
        "algo_block_size": str(block_size),
        "epsilon": "any",
        "similarity_func": "any",
        "scramble": "1"
    }
    
    for similarity in ["'scalar'","'jaccard'"]:
        fixed["similarity_func"] = similarity
        q = build_query(fixed)
        results_df.query(q).plot(x = "VBS_avg_nzblock_height", y = "relative_density", kind = "scatter", label = "similarity");
        
    plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
    #plt.xlim(0,1000)
    #plt.ylim(0,2)
    plt.axvline(block_size, alpha = 0.2, color = "red")
    plt.axhline(1, alpha = 0.2, color = "red")
    plt.title("M,N = {},{}\n block_size = {} \n density inside nonzero blocks = {} \n fraction of nonzero-blocks = {} \n similarity function = {}".format(cols,rows,block_size,  i_density, b_density, similarity))
    plt.xlabel("Average height of nonzero blocks");
    plt.ylabel("Average density inside nonzero blocks \n (relative to original blocking)");
    
    savename = save_folder + name + "_r{}_c{}_b{}_d{}_bd{}_{}.jpg".format(cols,rows,block_size,  i_density, b_density, similarity);

    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    plt.close()


def reorder_heatmap(rows = 2048, cols = 2048, block_size = 64, similarity = "'scalar'", save_folder = "../images/reorder_landscape/", name = "reorder_heatmap", saad = False):
    
    
    fixed = {
        "input_entries_density": None,
              "input_blocks_density": None,
              "rows": str(rows),
              "cols": str(cols),
              "input_block_size": str(block_size),
              "algo_block_size": str(block_size),
              "epsilon": "any",
              "similarity_func": similarity,
              "scramble": "1"
    }

    if saad: 
            fixed["algo_block_size"] = 1;
            name = "reorder_heatmap_saad";
    
    heatmap_array = []
    for input_blocks_density in results_df["input_blocks_density"].unique():
        
        row = []
        for input_entries_density in results_df["input_entries_density"].unique():
    
            fixed["input_entries_density"] = input_entries_density;
            fixed["input_blocks_density"] = input_blocks_density;
            q = build_query(fixed)
            interp_df = results_df.query(q).sort_values("VBS_avg_nzblock_height");
            interp_df.drop_duplicates("VBS_avg_nzblock_height", inplace = True)
            interp_df.sort_values("VBS_avg_nzblock_height", ascending = True, inplace = True)
            xp = interp_df.query(q)["VBS_avg_nzblock_height"];
            yp = interp_df.query(q)["relative_density"]
            interp = np.interp(64,xp,yp)
            row.append(interp)
        heatmap_array.append(row)
    
    
    heat_df = pd.DataFrame(heatmap_array, columns=results_df["input_entries_density"].unique())
    heat_df.set_index(results_df["input_blocks_density"].unique(), inplace = True)
    heat_df.sort_index(level=0, ascending=False, inplace=True)
    
    cmap = sns.diverging_palette(0,255,sep=1, as_cmap=True)
    
    plt.gca()
    ax = sns.heatmap(heat_df, linewidths=.5, annot=True, cbar_kws={'label': 'relative density'}, cmap = cmap, center = 1, vmin = 0, vmax = 1)
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)
    plt.xlabel("Density inside nonzero blocks") 
    plt.ylabel("Fraction of nonzero blocks");
    
    plt.title("M,N = {},{}\n block_size = {} \n similarity function = {}".format(cols,rows,block_size, similarity))
    savename = save_folder + name + "_r{}_c{}_b{}_{}.jpg".format(rows, cols, block_size, similarity);
    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    
    

for cols in [2048,]:
    for rows in [2048,]:
        for block_size in [16,32,64]:
            for similarity in ["'scalar'","'jaccard'"]:
                try:
                    reorder_heatmap(cols,rows,block_size, similarity, name = "reorder_heatmap", saad = saad);
                except:
                    print("could not make image for cols = {}, rows = {}, block_size = {}".format(cols,rows,block_size))