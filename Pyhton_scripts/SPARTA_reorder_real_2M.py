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

saving_path = "../images/reorder_experiment_real_2M";
try:
    os.mkdir(saving_path)
except:
    print("overwriting images");
    
datasetName = "../results/compression_results.txt"
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


def nonzero_blocks(i):
    return experiments["VBS_nz_blocks"][i]/(experiments["VBS_block_rows"][i]*experiments["cols"][i]/experiments["algo_block_size"][i]);

def VBS_density_inside_blocks_ratio(i):
    return (experiments["total_nonzeros"][i]/experiments['VBS_total_nonzeros'][i])/input_density(i);

def ncVBS_density_inside_blocks_ratio(i):
    return (experiments["total_nonzeros"][i]/ncVBS_total_nonzeros(i))/input_density(i)

def ncVBS_total_nonzeros(i):
    return experiments["ncVBS_height_mean"][i]*experiments["cols"][i];

def density_ratio(i):
    return ncVBS_density_inside_blocks_ratio(i)/VBS_density_inside_blocks_ratio(i);

def input_density(i):
    return experiments["total_nonzeros"][i]/(experiments["rows"][i]*experiments["cols"][i]);

def height_ratio(i):
    return experiments["VBS_block_rows"][i] /experiments['ncVBS_height_mean'][i]

def ncVBS_height(i):
    return experiments["rows"][i]/experiments["ncVBS_height_mean"][i]


def plot_VBS_density_vs_block_density(epsilon, block_size, scramble = 1):
    x_field_function = input_density;
    y_field_function = VBS_density_inside_blocks_ratio;
    fixed = {"algo_block_size": block_size, "scramble": scramble, "epsilon": epsilon};
    varying = "input_source";
    plot_x_vs_y_all_z(x_field_function, y_field_function, fixed = fixed, varying = varying);
    
    plt.title(make_title(fixed))
    plt.xlabel("input entries density");
    plt.ylabel("VBS in-block density multiplier");
    plt.xlim(0,0.005)
    plt.legend(bbox_to_anchor=(1, 1))
    
def plot_height_ratio(epsilon, block_size, scramble = 1):
    x_field_function = input_density;
    y_field_function = height_ratio;
    fixed = {"algo_block_size": block_size, "scramble": scramble, "epsilon": epsilon};
    varying = "input_source";
    plot_x_vs_y_all_z(x_field_function, y_field_function, fixed = fixed, varying = varying);
    
    plt.title(make_title(fixed))
    plt.xlabel("input entries density");
    plt.ylabel("ncVBS vs VBS blocks-height ratio");
    plt.xlim(0,0.005)
    plt.legend(bbox_to_anchor=(1, 1))

def plot_in_block_density_ratio(epsilon, block_size, scramble = 1):
    x_field_function = input_density;
    y_field_function = density_ratio;
    fixed = {"algo_block_size": block_size, "scramble": scramble, "epsilon": epsilon};
    varying = "input_source";
    plot_x_vs_y_all_z(x_field_function, y_field_function, fixed = fixed, varying = varying);
    
    plt.title(make_title(fixed))
    plt.xlabel("input entries density");
    plt.ylabel("ncVBS vs VBS in-block density ratio");
    plt.xlim(0,0.005)
    plt.legend(bbox_to_anchor=(1, 1))
    
    
def plot_ncVBS_density_multiplyer(epsilon, block_size, scramble = 1):
    x_field_function = input_density;
    y_field_function = ncVBS_density_inside_blocks_ratio;
    fixed = {"algo_block_size": block_size, "scramble": scramble, "epsilon": epsilon};
    varying = "input_source";
    plot_x_vs_y_all_z(x_field_function, y_field_function, fixed = fixed, varying = varying);
    
    plt.title(make_title(fixed))
    plt.xlabel("input entries density");
    plt.ylabel("ncVBS in-block density multiplier");
    plt.xlim(0,0.005)
    plt.legend(bbox_to_anchor=(1, 1))
    
        
scramble = 3;
epsilon = 0.7;
block_size = 128;

title_props = "s_" + str(scramble) + "_e_" + str(epsilon) + "_bsize_" + str(block_size) + ".jpg";

plt.figure()
plot_VBS_density_vs_block_density(epsilon = epsilon, block_size = block_size, scramble = scramble);
plt.savefig(saving_path + "/VBS_density_multiplier" + title_props);

plt.figure()
plot_height_ratio(epsilon = epsilon, block_size = block_size, scramble = scramble);
plt.savefig(saving_path + "/comparison_height" + title_props);

plt.figure()
plot_in_block_density_ratio(epsilon = epsilon, block_size = block_size, scramble = scramble);
plt.savefig(saving_path + "/comparison_density" + title_props);

plt.figure()
plot_ncVBS_density_multiplyer(epsilon = epsilon, block_size = block_size, scramble = scramble);
plt.savefig(saving_path + "/ncVBS_density_multiplier" + title_props);

