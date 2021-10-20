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

                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_reordering_blocked_synth_2.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="../images/",
        help="directory where the images are saved")

    args = parser.parse_args()

    input_csv = args.input_csv;
    output_dir = args.output_dir;
                    
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


to_display = ["input_entries_density",
              "input_blocks_density",
              "rows",
              "cols",
              "input_block_size",
              "algo_block_size",
              "epsilon",
              "similarity_func",
              "scramble",
        ]

for var in to_display:
    print(var, results_df[var].unique());



fixed = {
        "input_entries_density": "0.5",
              "input_blocks_density": "0.5",
              "rows": "2048",
              "cols": "2048",
              "input_block_size":"64",
              "algo_block_size":"64",
              "epsilon": "any",
              "similarity_func":"'scalar'",
              "scramble": "1"
}


def build_query(fixed):
    #build a query from the fixed dictionary
    q = ""
    for k, val in fixed.items():
        if val != "any":
            q += k + " == " + str(val) + " and ";
    q = q[0:-5];
    return q;


q = build_query(fixed)
results_df.query(q).size
results_df.query(q).plot(x = "VBS_avg_nzblock_height", y = "relative_density", kind = "scatter");
plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
plt.xlim(10,85)
plt.ylim(0.8,1.05)
plt.title("Reordering a scrambled blocked matrix (zoom in)")
#plt.axvline(64, alpha = 0.2, color = "red")
plt.xlabel("Average height of nonzero blocks");
plt.ylabel("Average density of nonzero blocks \n (relative to original blocking)");
plt.savefig(output_dir + "example_block_curve_detail.jpg", format = 'jpg', dpi=300, bbox_inches = "tight")



#TODO interpolate to find best density for given size. Produce heatmap


def reorder_heatmap(rows = 2048, cols = 2048, block_size = 64, similarity = "'scalar'", save_folder = "../images/reorders"):
    
    
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
    plt.ylabel("Density of nonzero blocks");
    
    plt.title("M,N = {},{}\n block_size = {} \n similarity function = {}".format(cols,rows,block_size, similarity))
    savename = save_folder + "reorder_heatmap_r{}_c{}_b{}_{}.jpg".format(rows, cols, block_size, similarity);
    plt.savefig(savename, format = 'jpg', dpi=300, bbox_inches = "tight")
    plt.show()
    
    

for cols in [2048,]:
    for rows in [2048,]:
        for block_size in [16,32,64]:
            for similarity in ["'scalar'", "'jaccard'"]:
                try:
                    reorder_heatmap(cols,rows,block_size, similarity);
                except:
                    print("could not make image for cols = {}, rows = {}, block_size = {}".format(cols,rows,block_size))




ax = plt.gca();
ax.set(ylabel="density of nonzero blocks", xlabel = "density inside blocks")
#q = "cols == 1024 and B_cols == 16384 and input_blocks_density == 64  and density < 0.1"
#results_df.query(q).plot(y = "input_blocks_density", x = "input_entries_density", kind = "scatter", c = "winner", colormap='viridis');


#ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlim([0.0005, 0.2]);
plt.savefig("landscape.jpg", format = 'jpg', dpi=300, bbox_inches = "tight")



plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "input_entries_density"

rows = 2048
cols = 2048;
B_cols = 8192;
input_block_size = 64;
input_blocks_density = 0.2

q = "rows == {} and cols == {} and B_cols == {} and input_block_size == {} and input_blocks_density == {}".format(rows, cols, B_cols, input_block_size, input_blocks_density)
results_df.query(q).plot(x = x, y = y, color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density inside blocks',
       title='Comparison with cusparse \n N, K, M = {}, {}, {} \n block size = {} \n input_blocks_density = {}'.format(rows, cols, B_cols, input_block_size, input_blocks_density))
ax.axhline(1)
x = results_df.query(q)["input_entries_density"].tolist();
y = np.array(results_df.query(q)["sp_vs_cu"].tolist());
errors = np.array(results_df.query(q)["VBSmm_std"].tolist());
upper = y + 3*errors;
lower = y - 3*errors;

ax.set_xscale("log")
ax.set_xlim([0.0005, 0.2]);

plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
         
ax.fill_between(x, lower, upper, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.savefig("sp_vs_cusparse_0.2.jpg", format = 'jpg', dpi=300, bbox_inches = "tight")


plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "input_entries_density"

rows = 2048
cols = 2048;
B_cols = 8192;
input_block_size = 64;
input_blocks_density = 0.5

q = "rows == {} and cols == {} and B_cols == {} and input_block_size == {} and input_blocks_density == {}".format(rows, cols, B_cols, input_block_size, input_blocks_density)
results_df.query(q).plot(x = x, y = y, color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density inside blocks',
       title='Comparison with cusparse \n N, K, M = {}, {}, {} \n block size = {} \n input_blocks_density = {}'.format(rows, cols, B_cols, input_block_size, input_blocks_density))
ax.axhline(1)
x = results_df.query(q)["input_entries_density"].tolist();
y = np.array(results_df.query(q)["sp_vs_cu"].tolist());
errors = np.array(results_df.query(q)["VBSmm_std"].tolist());
upper = y + 3*errors;
lower = y - 3*errors;

ax.set_xscale("log")
ax.set_xlim([0.0005, 0.2]);

plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
         
ax.fill_between(x, lower, upper, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.savefig("sp_vs_cusparse_0.5.jpg", format = 'jpg', dpi=300, bbox_inches = "tight")

