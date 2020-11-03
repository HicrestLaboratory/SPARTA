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
        if not check_good_sample(i):
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

def plot_x_vs_y_all_z(x_field_function, y_field_function, fixed, varying, draw_std_error = True):
    
    varying_values = find_unique(varying);
    for val in varying_values:
        new_fixed = dict(fixed);
        new_fixed[varying] = val;
        plot_x_vs_y(x_field_function, y_field_function, new_fixed, label = val);
    
    return 0

def make_title(fixed_values_dict):
    title = "";
    for key, value in fixed_values_dict.items():
        title += key + ":" + str(value) + "\n"
    return title

saving_path = "images/real_reorder_experiment";
try:
    os.mkdir(saving_path)
except:
    print("overwriting images");
    
datasetName = "real_reordering_results.txt"
experiments = {};



with open(datasetName, 'r') as f:
    fields = f.readline();
    fields = fields.split();
    for elem in fields:
        experiments[elem] = [];

    stop = False;
    while ( not stop ):
            line = f.readline();
            if not line:
                stop = True;
                break
            line = line.split();
            if (line[0] == "exp_name"):
                continue;
            for elem in zip(fields, line):
                data = elem[1];
                try:
                    data = float(data)
                except ValueError:
                    0;
                experiments[elem[0]].append(data);
            line = f.readline();
       
        
print(fields)      
            
print("IMPORT DONE")
            
n_exp = len(experiments["input_type"]);

def find_unique(property_name):
    value_list = np.array(experiments[property_name]);
    return np.unique(value_list);

graphs = find_unique("input_source");

def plot_density_vs_block_density(epsilon, block_size):
    plt.figure();
    x_field_function = lambda i: experiments["total_nonzeros"][i]/(experiments["rows"][i]*experiments["cols"][i]);
    y_field_function = lambda i: experiments["VBS_nz_blocks"][i]/(experiments["VBS_block_rows"][i]*experiments["cols"][i]/experiments["algo_block_size"][i]);
    fixed = {"algo_block_size": block_size, "scramble": 3, "epsilon": epsilon};
    varying = "input_source";
    plot_x_vs_y_all_z(x_field_function, y_field_function, fixed = fixed, varying = varying);
    
    plt.title(make_title(fixed))
    plt.xlabel("input entries density");
    plt.ylabel("percentage of nonzero blocks");
    plt.legend()
    plt.show()
    
    
plot_density_vs_block_density(epsilon = 0.8, block_size = 32);