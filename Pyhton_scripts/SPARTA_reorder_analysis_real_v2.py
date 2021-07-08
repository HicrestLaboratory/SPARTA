# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:04:46 2021

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
name = "reordering_vs_saad_fixed"
filename = folder + name + ".txt";
csv_filename = folder + name + ".csv";
                    
convert_to_CSV(filename, csv_filename);
results_df = pd.read_csv(csv_filename);

data_top = results_df.columns 

columsn = ['exp_name', 
           'input_type', 
           'input_source', 
           'reorder_algorithm', 
           'rows',
           'cols',
           'total_nonzeros',
           'input_blocks_density',
           'input_entries_density',
           'input_block_size',
           'algo_block_size',
           'epsilon',
           'scramble',
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

results_df["density"] = results_df.apply(lambda x: x['total_nonzeros']/(x['rows']*x["cols"]), axis=1)
results_df["VBS_density"] = results_df.apply(lambda x: x['total_nonzeros']/x['VBS_total_nonzeros'], axis=1)


source = "data/real_graphs/CA-CondMat.txt";

plt.figure()

scramble = 1;
algo_block_size = 32;
rows = 6110;

ax = plt.gca();
x = "VBS_avg_nzblock_height";
y = "VBS_density"
q = 'rows == @rows and algo_block_size == @algo_block_size and scramble == @scramble and reorder_algorithm == "saad"' 
results_df.query(q).plot(x = x, y = y, kind = "scatter", ax = ax);


plt.title("Saad's blocking for facebook SNAP graph \n block width fixed at 32")
plt.ylabel("density after blocking");
plt.xlabel("average block height");

#q = 'rows == @rows and algo_block_size == @algo_block_size and scramble == @scramble and reorder_algorithm == "saad_blocks"' 
#results_df.query(q).plot(x = x, y = y, kind = "scatter", color = "red", ax = ax);
