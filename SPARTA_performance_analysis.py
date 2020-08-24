# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 12:42:17 2020

@author: Paolo
"""
import os
import matplotlib.pyplot as plt
import numpy as np


saving_path = "images/reorder_experiment";
try:
    os.mkdir(saving_path)
except:
    print("overwriting images");
    
datasetName = "test_reordering.txt"
experiments = {};

with open(datasetName, 'r') as f:
    fields = f.readline();
    fields = fields.split();
    for elem in fields:
        experiments[elem] = [];

    while (True):
        line = f.readline();
        if not line:
            break;
        line = line.split();
        if (line[0] == "input_type"):
            continue;
        for elem in zip(fields, line):
            experiments[elem[0]].append(float(elem[1]));
            
            
            
n_exp = len(experiments["input_type"]);


def npslice(var_name,eps,algo_bs):
    return np.array(experiments[var_name])[i_vals[eps][algo_bs]];

input_block_size = 64;
input_block_density = 0.5;
i_vals = {};


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
        plt.xlabel("density inside input blocks");
        plt.ylabel("output fill-in ratio \n (higher is better)");
        plt.title("epsilon = " + str(eps) + "\n" + "input block size = " + str(input_block_size) + \
                  "\n input density of nonzero blocks = " + str(input_block_density)\
                  );
        axes.set_ylim(0,1);
        axes.set_xlim(0,1);

    plt.tight_layout()
    plt.legend(title = "algorithm block size");
    plt.savefig(saving_path + '/density_vs_fillin_varying_blocksize_eps' + str(eps) + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
    plt.savefig(saving_path + '/density_vs_fillin_varying_blocksize_eps' + str(eps) + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")


plt.figure();
algo_bs = 50;
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
plt.xlabel("density inside blocks");
plt.ylabel("fill-in ratio \n (higher is better)");
plt.title(title)
plt.legend(title = "epsilon");
plt.savefig(saving_path + '/density_vs_fillin_varying_epsilon.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/density_vs_fillin_varying_epsilon.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")




input_block_size = 64;
input_entries_density = 0.8 ;
ibd = [0.1, 0.2, 0.4]
algo_block_size = 35;

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
plt.ylabel("vertex compression ratio \n (higher is better)");
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
    plt.ylabel("output fill in ratio \n (higher is better)");
    plt.plot(x,y, label = input_block_density);
plt.tight_layout()
plt.legend(title = "nonzero block density");
plt.title(title);
plt.savefig(saving_path + '/epsilon_vs_fillin' + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
plt.savefig(saving_path + '/epsilon_vs_fillin' + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")




