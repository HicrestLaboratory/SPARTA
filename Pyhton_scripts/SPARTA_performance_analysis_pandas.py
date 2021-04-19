# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 11:16:49 2021

@author: Paolo
"""


import pandas as pd
import os
import matplotlib.pyplot as plt

def change_sep(string, newsep, oldsep = " "):
    string = string.split(oldsep);
    newstring = "";
    for elem in string:
        newstring += elem + newsep
    return newstring[:-1];

def convert_to_CSV(filename, new_file_name):
    with open(filename, 'r') as infile:
        with open(new_file_name, 'w') as outfile:
            fieldsline = infile.readline();
            
            firstrow = change_sep(fieldsline, ",");
            
            outfile.writelines(firstrow);
        
            stop = False;
            while ( not stop ):
                    line = infile.readline();
                    if not line:
                        stop = True;
                        break
                    newline = change_sep(line, ",")
                    outfile.writelines(newline);                    
                    line = infile.readline();
                    
                    
folder = "../results/";
name = "test_complete_cublas_results_ultrasparse_5"
filename = folder + name + ".txt";
csv_filename = folder + name + ".csv";
                    
convert_to_CSV(filename, csv_filename);
results_df = pd.read_csv(csv_filename);

columns = ['input_type',
 'A_rows',
 'A_cols',
 'A_total_nonzeros',
 'A_blocks_density',
 'A_entries_density',
 'A_block_size',
 'B_cols',
 'B_density',
 'Warmup',
 'Repetitions',
 'Algorithm',
 'gemm_mean(ms)',
 'gemm_std',
 'VBSmm_mean(ms)',
 'VBSmm_std',
 'cusparse_spmm_mean(ms)',
 'cusparse_spmm_std']


results_df["density"] = results_df.apply(lambda x: x['A_total_nonzeros']/(x['A_rows']*x["A_cols"]), axis=1)
results_df["sp_vs_cu"] = results_df.apply(lambda x: x['cusparse_spmm_mean(ms)']/x["VBSmm_mean(ms)"], axis=1)
results_df["sp_vs_gemm"] = results_df.apply(lambda x: x['gemm_mean(ms)']/x["VBSmm_mean(ms)"], axis=1)

results_df["in-block-density"] = results_df.apply(lambda x: x['density']/x["A_blocks_density"], axis=1)
def winner(x):
    if (x["VBSmm_mean(ms)"] < x["cusparse_spmm_mean(ms)"] and x["VBSmm_mean(ms)"] < x["gemm_mean(ms)"]):
        return 0;
    elif (x["gemm_mean(ms)"] < x["cusparse_spmm_mean(ms)"]):
        return 1;
    else:
        return 2;

results_df["winner"] = results_df.apply(winner, axis = 1);
q = "A_cols == 1024 and B_cols == 16384 and A_block_size == 128 and density < 0.1"
results_df.query(q).plot(y = "A_blocks_density", x = "A_entries_density", kind = "scatter", c = "winner", colormap='viridis');

q = "A_cols == 1024 and B_cols == 16384 and A_block_size == 64  and density < 0.1"
results_df.query(q).plot(y = "A_blocks_density", x = "A_entries_density", kind = "scatter", c = "winner", colormap='viridis');


plt.figure()
ax = plt.gca();
y = "sp_vs_gemm";
x = "A_blocks_density"
q = "A_cols == 1024 and B_cols == 16384 and A_block_size == 256"
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "green", ax = ax);


q = "A_cols == 2048 and B_cols == 16384 and A_block_size == 128"
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);


q = "A_cols == 2048 and B_cols == 16384 and A_block_size == 64"
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "blue", ax = ax);



plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "A_blocks_density"

A_cols = 1024;
B_cols = 16384;
A_block_size = 128;
density = 0.005

q = "A_cols == {} and B_cols == {} and A_block_size == {} and density < {}".format(A_cols, B_cols, A_block_size, density)
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density of blocks',
       title='Comparison with cusparse \n N, K, M = 1024, {}, {} \n block size = {} \n global density = {}'.format(A_cols, B_cols, A_block_size, density))
ax.axhline(1)

plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "B_cols"

A_cols = 1024;
A_blocks_density = 0.15;
A_block_size = 128;
density = 0.005

q = "A_cols == {} and A_block_size == {} and A_blocks_density == {} and density < {}".format(A_cols, A_block_size, A_blocks_density, density)
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='M',
       title='Comparison with cusparse \n N, K,= 1024, {}, \n density of blocks = {} \n block size = {} \n global density = {}'.format(A_cols, A_blocks_density, A_block_size, density))
ax.axhline(1)

plt.figure()
ax = plt.gca();
y = "sp_vs_cu";
x = "A_blocks_density"

A_cols = 1024;
B_cols = 16384;
A_block_size = 64;
density = 0.001

q = "A_cols == {} and B_cols == {} and A_block_size == {} and density < {}".format(A_cols, B_cols, A_block_size, density)
results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
ax.set(ylabel='speedup vs cusparse', xlabel='density of blocks',
       title='Comparison with cusparse \n N, K, M = 1024, {}, {} \n block size = {} \n global density = {}'.format(A_cols, B_cols, A_block_size, density))
ax.axhline(1)


