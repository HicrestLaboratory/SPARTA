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
    
    parser.add_argument("--input-csv", default="../results/test_cublas_results_ultrasparse.csv",
        help="file that contains the already run experiments")
    parser.add_argument("--output-dir", default="",
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
    
print(results_df["input_entries_density"].unique())
print(results_df["input_blocks_density"].unique())


results_df["winner"] = results_df.apply(winner, axis = 1);
q = "cols == 1024 and rows == 1024 and B_cols == 4096 and input_block_size == 64"
results_df.query(q).plot(y = "input_blocks_density", x = "input_entries_density", kind = "scatter", c = "winner", colormap='viridis');


plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)

         
         
ax = plt.gca();
#q = "cols == 1024 and B_cols == 16384 and input_blocks_density == 64  and density < 0.1"
#results_df.query(q).plot(y = "input_blocks_density", x = "input_entries_density", kind = "scatter", c = "winner", colormap='viridis');


#ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlim([0.0005, 0.3]);
plt.show();


plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "input_entries_density"

rows = 2048
cols = 4096;
B_cols = 4096;
input_block_size = 64;
input_blocks_density = 0.5

q = "rows == {} and cols == {} and B_cols == {} and input_block_size == {} and input_blocks_density == {}".format(rows, cols, B_cols, input_block_size, input_blocks_density)
results_df.query(q).plot(x = x, y = y, color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density of blocks',
       title='Comparison with cusparse \n N, K, M = {}, {}, {} \n block size = {} \n input_blocks_density = {}'.format(rows, cols, B_cols, input_block_size, input_blocks_density))
ax.axhline(1)
x = results_df.query(q)["input_entries_density"].tolist();
y = np.array(results_df.query(q)["sp_vs_cu"].tolist());
errors = np.array(results_df.query(q)["VBSmm_std"].tolist());
upper = y + 3*errors;
lower = y - 3*errors;

plt.grid(b=True, which='major', color='#666666', linestyle='-', alpha=0.4)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2, linewidth=0.5)
         
ax.fill_between(x, lower, upper, alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.show()



plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "B_cols"

cols = 1024;
input_blocks_density = 0.15;
input_blocks_density = 128;
density = 0.5

q = "cols == {} and input_blocks_density == {} and input_blocks_density == {} and density < {}".format(cols, input_blocks_density, input_blocks_density, density)
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='M',
       title='Comparison with cusparse \n N, K,= 1024, {}, \n density of blocks = {} \n block size = {} \n global density = {}'.format(cols, input_blocks_density, input_blocks_density, density))
ax.axhline(1)

plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "input_blocks_density"

cols = 1024;
B_cols = 16384;
input_blocks_density = 64;
density = 0.001

q = "cols == {} and B_cols == {} and input_blocks_density == {} and density < {}".format(cols, B_cols, input_blocks_density, density)
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density of blocks',
       title='Comparison with cusparse \n N, K, M = 1024, {}, {} \n block size = {} \n global density = {}'.format(cols, B_cols, input_blocks_density, density))
ax.axhline(1)


