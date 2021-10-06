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
                    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Makes images from experiments")
    
    parser.add_argument("--input-csv", default="../results/test_reordering_blocked_synth.csv",
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
        "input_entries_density": "0.1",
              "input_blocks_density": "0.1",
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

results_df.query(q).plot(x = "relative_density", y = "output_in_block_density", kind = "scatter");
plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
#plt.xlim(0,64)
plt.savefig(output_dir + "block_curve.jpg", format = 'jpg', dpi=300, bbox_inches = "tight")


         
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

