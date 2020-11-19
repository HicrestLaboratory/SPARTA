# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 12:42:17 2020

@author: Paolo
"""
import os
import matplotlib.pyplot as plt
import numpy as np

def check_fixed_condition(val, condition):
    if type(condition) is list:
        return (val >= condition[0] and val <= condition[1]);
    else:
        return val == condition;
    
def check_good_sample(i):
    return (experiments["total_nonzeros"][i] != 0);


def plot_x_vs_y(x_field_function, y_field_function, fixed, label = None, draw_std_error = True):
    y_lists = {};
    y_values = [];
    y_errors =[];
    for i in range(n_exp):
        #check condition for adding experiment
        skip = False; 
        for fixed_field, condition in fixed.items():
            if not check_fixed_condition(experiments[fixed_field][i], condition):
                skip = True;
            if not check_good_sample:
                skip = True;
        if skip:
            continue;   
        
        x = x_field_function(i);
        if x not in y_lists:
            y_lists[x] = [];
        
        y = y_field_function(i);
        
        #append to data list
        y_lists[x].append(y);
        
    for x in y_lists.keys():
        y_values.append(np.mean(y_lists[x]));
        y_errors.append(np.std(y_lists[x]));

    x_values = list(y_lists.keys());
    plt.scatter(x_values,y_values, label = label);
    
    if draw_std_error:
        plt.fill_between(x_values, np.array(y_values) - np.array(y_errors),
                             np.array(y_values) + np.array(y_errors), alpha=0.1,
                             color="r");



saving_path = "images/reorder_experiment_v2";
try:
    os.mkdir(saving_path)
except:
    print("overwriting images");
    
datasetName = "test_reordering_old.txt"
experiments = {};

with open(datasetName, 'r') as f:
    fields = f.readline();
    fields = fields.split();
    for elem in fields:
        experiments[elem] = [];
    experiments["scrambled"] = [];

    stop = False;
    while ( not stop ):
        for scrambled in [0,1]:
            line = f.readline();
            if not line:
                stop = True;
                break
            line = line.split();
            if (line[0] == "input_type"):
                continue;
            if (scrambled == 0):
                for elem in zip(fields, line):
                    experiments[elem[0]].append(float(elem[1]));
                experiments["scrambled"].append(scrambled); #s is missed in the output report, so added here
            line = f.readline();
       
        
print(fields)
def npslice(var_name,eps,algo_bs):
    return np.array(experiments[var_name])[i_vals[eps][algo_bs]];            
            
print("IMPORT DONE")
            
n_exp = len(experiments["input_type"]);

input_block_size = 64;
input_block_density = 0.01;
i_vals = {};


plt.figure();
x_field_function = lambda i: experiments["input_entries_density"][i]*experiments["input_blocks_density"][i];
y_field_function = lambda i: experiments["total_nonzeros"][i]/experiments["VBS_total_nonzeros"][i];
for algo_block_size in [100,]:
    fixed = {"input_blocks_density": 0.9, "input_block_size": 64, "algo_block_size": algo_block_size, "epsilon" : 0.99};
    plot_x_vs_y(x_field_function, y_field_function, fixed, label = algo_block_size);

plt.xlabel("input density");
plt.ylabel("output in-block density");
plt.legend()
plt.show()



plt.figure()
for i in range(n_exp):
    if experiments["total_nonzeros"][i] != 0 and\
    experiments["input_blocks_density"][i] == input_block_density and\
    experiments["input_block_size"][i] == input_block_size:

        
        ev = experiments["epsilon"][i];
        if ev not in i_vals:
            i_vals[ev] = {};
        
        algo_bs = experiments["algo_block_size"][i];
        if algo_bs not in i_vals[ev]:
            i_vals[ev][algo_bs] = [];

        i_vals[ev][algo_bs].append(i);

for eps in i_vals:
    
    plt.figure();
    for algo_bs in i_vals[eps]:
        efficiency = [];
        for i in i_vals[eps][algo_bs]:
            nz = experiments["total_nonzeros"][i];
            vbs_nz = experiments["VBS_total_nonzeros"][i];
            efficiency.append(nz/vbs_nz);
    
        plt.scatter(npslice("input_entries_density",eps,algo_bs),efficiency, label = str(algo_bs));
        axes = plt.gca();
        plt.xlabel("input in-block density");
        plt.ylabel("output in-block density \n (higher is better)");
        plt.title("epsilon = " + str(eps) + "\n" + "input block size = " + str(input_block_size) + \
                  "\n input density of nonzero blocks = " + str(input_block_density)\
                  );
        axes.set_ylim(0,1);
        axes.set_xlim(0,1);

    plt.tight_layout()
    plt.legend(title = "algorithm block size");
    plt.savefig(saving_path + '/Idensity_vs_Odensity_varying_blocksize_eps' + str(eps) + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
    plt.savefig(saving_path + '/Idensity_vs_Odensity_varying_blocksize_eps' + str(eps) + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")


plt.figure();
algo_bs = 11;
title = "algorithm block size = " + str(algo_bs) +\
    "\n input block size = " + str(input_block_size) +\
    "\n input density of nonzero blocks = " + str(input_block_density);

for eps in i_vals:
    efficiency = [];
    for i in i_vals[eps][algo_bs]:
        nz = experiments["total_nonzeros"][i];
        vbs_nz = experiments["VBS_total_nonzeros"][i];
        efficiency.append(nz/vbs_nz);
    plt.scatter(npslice("input_entries_density",eps,algo_bs),efficiency, label = str(eps));
plt.tight_layout()
plt.xlabel("input in-block density");
plt.ylabel("output in-block density \n (higher is better)");
plt.title(title)
plt.legend(title = "epsilon");
plt.savefig(saving_path + '/Idensity_vs_Odensity_varying_epsilon.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/Idensity_vs_Odensity_varying_epsilon.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")




input_block_size = 64;
input_entries_density = 0.2 ;
ibd = [0.01, 0.2, 0.4]
algo_block_size = 37;

#    "\n input density of nonzero blocks = " + str(input_block_density) +\

title = "algorithm block size = " + str(algo_block_size) +\
    "\n input block size = " + str(input_block_size) +\
    "\n input density inside blocks = " + str(input_entries_density);

plt.figure();
for input_block_density in ibd:
    exp2 = {};
    for i in range(n_exp):
        if experiments["total_nonzeros"][i] != 0 and\
        experiments["input_blocks_density"][i] == input_block_density and\
        experiments["input_block_size"][i] == input_block_size and\
        experiments["input_entries_density"][i] == input_entries_density and\
        experiments["algo_block_size"][i] == algo_block_size:
            epsilon = experiments["epsilon"][i];
            if epsilon not in exp2:
                exp2[epsilon] = [];
            rows = experiments["rows"][i];
            block_rows = experiments["VBS_block_rows"][i];
            compression_rate = rows/block_rows;
            exp2[epsilon].append(compression_rate);
    
    for eps in exp2:
        exp2[eps] = np.mean(exp2[eps]);
    
    lists = sorted(exp2.items());
    x,y = zip(*lists);    
    
    plt.plot(x,y, label = input_block_density);
 
plt.tight_layout()
plt.xlabel("epsilon");
plt.ylabel("average size of output block-rows \n (higher is better)");
plt.title(title);
plt.legend(title = "nonzero block density")
plt.savefig(saving_path + '/epsilon_vs_comp_ratio' + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/epsilon_vs_comp_ratio' + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")




plt.figure();
for input_block_density in ibd:
    exp2 = {};
    for i in range(n_exp):
        if experiments["total_nonzeros"][i] != 0 and\
        experiments["input_blocks_density"][i] == input_block_density and\
        experiments["input_block_size"][i] == input_block_size and\
        experiments["input_entries_density"][i] == input_entries_density and\
        experiments["algo_block_size"][i] == algo_block_size:
            epsilon = experiments["epsilon"][i];
            if epsilon not in exp2:
                exp2[epsilon] = [];
            rows = experiments["rows"][i];
            nz = experiments["total_nonzeros"][i];
            algo_nz = experiments["VBS_total_nonzeros"][i];
    
            fill_in = (nz/algo_nz);
            exp2[epsilon].append(fill_in);
    
    for eps in exp2:
        exp2[eps] = np.mean(exp2[eps]);
    
    lists = sorted(exp2.items());
    x,y = zip(*lists);    
    
    plt.xlabel("epsilon");
    plt.ylabel("output in-block density \n (higher is better)");
    plt.plot(x,y, label = input_block_density);
plt.tight_layout()
plt.legend(title = "nonzero block density");
plt.title(title);
plt.savefig(saving_path + '/epsilon_vs_Odensity' + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/epsilon_vs_Odensity' + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")




