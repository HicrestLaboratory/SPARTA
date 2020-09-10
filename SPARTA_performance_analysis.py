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


def add_column(field_name, function):
    experiments[field_name] = [];
    for i in range(n_exp):
        experiments[field_name].append(function(i));

def sum_std_err(errors):
    result = 0;
    N = len(errors);
    for s in errors:
        result += s**2;
    result = np.sqrt(result);
    return result/N;    


def clean_data(check):
    idx_list = [];    
    for i in range(n_exp):
        if check(i):
            idx_list.append(i);
    return idx_list;

def is_valid(i):
    if experiments["A_total_nonzeros"][i] != 0.:
        return True;
    
def mat_size(i):
    return experiments["A_cols"][i]*experiments["A_rows"][i];

def mat_sparsity(i):
    return experiments["A_total_nonzeros"][i]/mat_size(i)
    
def plot_best_scatter(x_field_func, y_field_func, method, fixed):
    point_list = {};
    x_list = [];
    y_list = [];
    val_list = [];
    methods = ["VBSmm_mean(ms)", "gemm_mean(ms)", "cusparse_spmm_mean(ms)"];
    for i in range(n_exp):
        skip = False; 
        for fixed_field, condition in fixed.items():
            if not check_fixed_condition(experiments[fixed_field][i], condition):
                skip = True;
        if skip:
            continue;
        x = x_field_function(i);
        y = y_field_function(i);
        if (x,y) not in point_list:
            point_list[(x,y)] = []

        min_met = min([experiments[met][i] for met in methods]);
        if min_met == experiments[method][i]:
            point_list[(x,y)].append(1);
        else:
            point_list[(x,y)].append(0);
            
    for k1,k2 in point_list.keys():
        avg = np.mean(point_list[(k1,k2)])
        if avg > 0.5:
            x_list.append(k1);
            y_list.append(k2);
            val_list.append(avg);
    plt.scatter(x_list, y_list, label = method);
    print(len(x_list));
    return 1

def plot_times_vs_other(x_field_function, time_field_function, time_std_field, fixed, label = None, draw_std_error = False):
    y_lists = {};
    y_lists_std = {};
    y_values = [];
    y_errors = [];
    for i in range(n_exp):
        #check condition for adding experiment
        skip = False; 
        for fixed_field, condition in fixed.items():
            if not check_fixed_condition(experiments[fixed_field][i], condition):
                skip = True;
        if skip:
            continue;
        
        A_cols = experiments["A_cols"][i];
        B_cols = experiments["B_cols"][i];
        A_total_nonzeros = experiments["A_total_nonzeros"][i];        
        
        x = x_field_function(i);
        if x not in y_lists:
            y_lists[x] = [];
            y_lists_std[x] = [];
        
        #convert time and its error to operations/ms
        total_multiplications = A_total_nonzeros*A_cols*B_cols;
        
        #y = total_multiplications/experiments[time_field][i];
        y = time_field_function(i);
        #y_std = (total_multiplications**2)/experiments[time_std_field][i];
        y_std = experiments[time_std_field][i];
        
        #append to data list
        y_lists[x].append(y);
        y_lists_std[x].append(y_std);
        
    for x in y_lists.keys():
        y_values.append(np.mean(y_lists[x]));
        y_err = sum_std_err(y_lists_std[x]);
        y_errors.append(y_err);

    x_values = list(y_lists.keys());
    print(len(x_values), len(y_values));
    print(x_values[:10])
    plt.scatter(x_values,y_values, label = label);
    
    if draw_std_error:
        plt.fill_between(x_values, np.array(y_values) - np.array(y_errors),
                             np.array(y_values) + np.array(y_errors), alpha=0.1,
                             color="r");



 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 #*****************************************************************************************************
 
 
 
save = True;
if save:
    saving_path = "images/performance_experiment_v2";
    try:
        os.mkdir(saving_path)
    except:
        print("overwriting images");
    
    
field_names = "input_type A_rows A_cols A_total_nonzeros A_blocks_density A_block_size \
                B_cols B_density VBS_AAM_nonzeros VBS_AAM_block_rows VBS_AAM_nz_blocks VBS_AAM_min_block_H VBS_AAM_max_block_H \
                Warmup Repetitions Algorithm gemm_mean(ms) gemm_std \
                VBSmm_mean(ms) VBSmm_std \
                VBSmm_nozeros_mean(ms) VBSmm_nozeros_std \
                VBSmm_angle_mean(ms) \
                VBSmm_angle_std cusparse_spmm_mean(ms) cusparse_spmm_std"    
    
datasetName = "test_complete_cublas_results.txt"
experiments = {};

expanded_names = {"A_rows" : "M", "A_cols": "K", "B_cols": "N", \
                  "A_total_nonzeros" : "total nonzeros in A", "A_blocks_density" : "percentage of nonzero blocks", \
                  "A_block_size": "block size", "A_in_block_density" : "in-block density"};


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
clean_idxs = clean_data(is_valid);
for key in experiments.keys():
    experiments[key] = np.array(experiments[key]);
    experiments[key] = experiments[key][clean_idxs];
    
n_exp = len(experiments["input_type"]);
print("num exp", n_exp);

#ADDING NEW FIELDS
#----------------------------------------------------------------------------------
def entries_density(i):
    A_rows = experiments["A_rows"][i];
    A_cols = experiments["A_cols"][i];
    A_total_nonzeros = experiments["A_total_nonzeros"][i];
    res = ((A_total_nonzeros)/(A_rows*A_cols));
    return  res;
add_column("A_entries_density", entries_density);

def in_block_density(i):
    A_rows = experiments["A_rows"][i];
    A_cols = experiments["A_cols"][i];
    A_total_nonzeros = experiments["A_total_nonzeros"][i];
    A_blocks_density = experiments["A_blocks_density"][i];
    res = (1/A_blocks_density)*((A_total_nonzeros)/(A_rows*A_cols));
    return  res;
add_column("A_in_block_density", in_block_density);

def relative_block_size(i):
    A_cols = experiments["A_rows"][i];
    relative_block_size = experiments["A_block_size"][i]/A_cols;
    return relative_block_size;
add_column("A_relative_block_size", relative_block_size);
#----------------------------------------------------------------------------------
 

dothis = False;
if dothis:
    x_field = "A_in_block_density";
    algo_name = "cusparse_spmm" 
    y_field = algo_name + "_mean(ms)";
    y_std_field = algo_name + "_std";
    
    x_field_function = lambda i: experiments[x_field][i]
    #x_field_function = relative_block_size;
    
    y_field_function = lambda i: experiments["cusparse_spmm_mean(ms)"][i]/experiments["VBSmm_mean(ms)"][i];
    
    
    #fixed = {"A_cols"  : 2048, "A_rows" : 2048, "B_cols" : 1024, "A_block_size" : 128, "A_in_block_density" : [0.5,0.6]};
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096*4; 
    block_density = 0.5;
    plot_name = "VBS_vs_spmm_A"+ str(A_cols) + "_B_" + str(B_cols) + "_fixed_blockdensity_" + str(block_density)
    title = "percentage of nonzero blocks = " + str(block_density) + "\n" + "M,K,N = (" + str(A_rows) + "," + str(A_cols) + "," + str(B_cols) + ")"; 
    
    for block_size in [128,256,512]:
        fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_blocks_density": block_density, "A_block_size" : block_size};
        plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = block_size);    
    
    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("VBS speedup vs cusparse \n (higher is better)");
    
    #plt.ylim([0,2]);
    #plt.xlim([0,0.3]);
    plt.legend(title = "Block size")
    plt.title(title);
    if save: 
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")





#--------------------------------------------------------------------------------------------------------

dothis = False;
if dothis:
    x_field = "A_in_block_density";
    algo_name = "cusparse_spmm" 
    y_field = algo_name + "_mean(ms)";
    y_std_field = algo_name + "_std";
    
    x_field_function = lambda i: experiments[x_field][i]
    #x_field_function = relative_block_size;
    
    y_field_function = lambda i: experiments["cusparse_spmm_mean(ms)"][i]/experiments["VBSmm_mean(ms)"][i];
    
    
    #fixed = {"A_cols"  : 2048, "A_rows" : 2048, "B_cols" : 1024, "A_block_size" : 128, "A_in_block_density" : [0.5,0.6]};
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096*4; 
    block_density = 0.3;
    block_size = 128
    plot_name = "VBS_vs_spmm_A"+ str(A_cols) + "_Bdensity_" + str(block_density) + "Bsize_" + str(block_size) + "_varying_MATsize";
    title = "percentage of nonzero blocks = " + str(block_density) + "\n" + "M,K = (" + str(A_rows) + "," + str(A_cols) + ")" + " \n Block size = " + str(block_size);
    
    for B_cols in [1024,2048,4096,4096*2, 4096*4]:
        fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_blocks_density": block_density, "A_block_size" : block_size};
        plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = B_cols);    
    
    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("VBS speedup vs cusparse \n (higher is better)");
    
    #plt.ylim([0,2]);
    #plt.xlim([0,1]);
    plt.legend(title = "Columns in B")
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")



#VBS vs cusparse, in_block_density. 
dothis = False
if dothis:
    x_field = "A_in_block_density";
    algo_name = "cusparse_spmm" 
    y_field = algo_name + "_mean(ms)";
    y_std_field = algo_name + "_std";
    
    x_field_function = lambda i: experiments[x_field][i]
    #x_field_function = relative_block_size;
    
    y_field_function = lambda i: experiments["cusparse_spmm_mean(ms)"][i]/experiments["VBSmm_mean(ms)"][i];
    
    
    #fixed = {"A_cols"  : 2048, "A_rows" : 2048, "B_cols" : 1024, "A_block_size" : 128, "A_in_block_density" : [0.5,0.6]};
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096; 
    block_density = 0.3;
    block_size = 128
    plot_name = "VBS_vs_spmm_A"+ str(A_cols) + "_B_" + str(B_cols) + "Block_size_" + str(block_size) + "_varying_Block_density";
    title = "\n" + "M,K,N = (" + str(A_rows) + "," + str(A_cols) + "," + str(B_cols) + ")" + " \n Block size = " + str(block_size);
    
    for block_density in [0.1,0.5,0.8]:
        fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_blocks_density": block_density, "A_block_size" : block_size};
        plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = block_density);    
    
    
    plt.axhline(y=1.0, color='r', linestyle='-', alpha = 0.5)

    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("VBS speedup vs cusparse \n (higher is better)");
    
    #plt.ylim([0,1.5]);
    #plt.xlim([0,0.3]);
    plt.legend(title = "Percentage of nonzero blocks")
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")

#VBS vs cusparse, in_block_density. ZOOMED 
dothis = False
if dothis:
    x_field = "A_in_block_density";
    algo_name = "cusparse_spmm" 
    y_field = algo_name + "_mean(ms)";
    y_std_field = algo_name + "_std";
    
    x_field_function = lambda i: experiments[x_field][i]
    #x_field_function = relative_block_size;
    
    y_field_function = lambda i: experiments["cusparse_spmm_mean(ms)"][i]/experiments["VBSmm_mean(ms)"][i];
    
    
    #fixed = {"A_cols"  : 2048, "A_rows" : 2048, "B_cols" : 1024, "A_block_size" : 128, "A_in_block_density" : [0.5,0.6]};
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096; 
    block_density = 0.3;
    block_size = 128
    plot_name = "VBS_vs_spmm_A"+ str(A_cols) + "_B_" + str(B_cols) + "Block_size_" + str(block_size) + "_varying_Block_density_ZOOMED";
    title = "\n" + "M,K,N = (" + str(A_rows) + "," + str(A_cols) + "," + str(B_cols) + ")" + " \n Block size = " + str(block_size);
    
    for block_density in [0.1, 0.5,0.8]:
        fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_blocks_density": block_density, "A_block_size" : block_size};
        plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = block_density);    


    plt.xlim([0,0.26]);
    plt.ylim([0,8]);
    
    plt.axhline(y=1.0, color='r', linestyle='-',alpha = 0.5)
    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("VBS speedup vs cusparse \n (higher is better)");
    
    plt.legend(title = "Percentage of nonzero blocks")
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")



#ABSOLUTE VBS vs cusparse, varying in-block density, fixed blocksize and block density
dothis = False
if dothis:
    x_field = "A_in_block_density";
    algo_name = "cusparse_spmm" 
    y_field = algo_name + "_mean(ms)";
    y_std_field = algo_name + "_std";
    
    x_field_function = lambda i: experiments[x_field][i] 

    #x_field_function = relative_block_size;    
    
    #fixed = {"A_cols"  : 2048, "A_rows" : 2048, "B_cols" : 1024, "A_block_size" : 128, "A_in_block_density" : [0.5,0.6]};
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096*4; 
    block_density = 0.2
    block_size = 512
    plot_name = "VBS_vs_spmm_absolute_A"+ str(A_cols) + "_B_" + str(B_cols) + "Block_size_" + str(block_size) + "_varying_in_block_density";
    title = "\n" + "M,K,N = (" + str(A_rows) + "," + str(A_cols) + "," + str(B_cols) + ")" + " \n Block size = " + str(block_size) + " \n percentage of nonzero blocks: " + str(block_density);
    
    
    y_fields = ["VBSmm_mean(ms)","gemm_mean(ms)"];
        
    y_field_function = lambda i: experiments[y_field][i]

    fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_blocks_density": block_density, "A_block_size" : block_size};
    plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = y_field);    
    
    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("time (ms) \n (smaller is better)");
    

    plt.legend(title = "Percentage of nonzero blocks")
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")



#VBS vs GEMM varying columns
dothis = False
if dothis:
    x_field = "A_blocks_density";
    x_field_function = lambda i: experiments[x_field][i];

    y_field_function = lambda i: experiments["gemm_mean(ms)"][i]/experiments["VBSmm_mean(ms)"][i];
    
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096*4; 
    in_block_density = 0.7;
    block_size = 128
    plot_name = "VBS_vs_gemm_A"+ str(A_cols) + "_B_" + str(B_cols) + "Block_size_" + str(block_size) + "_varying_Block_density";
    title = "\n" + "M,K = (" + str(A_rows) + "," + str(A_cols) + ")" + " \n Block size = " + str(block_size);
    
    
    for B_cols in [4096, 4096*2, 4096*4]:
        fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_block_size" : block_size};
        plot_times_vs_other(x_field_function, y_field_function, y_std_field, fixed, label = "N =" + str(B_cols));    

    plt.axhline(y=1.0, color='r', linestyle='-', alpha = 0.5)
    plt.tight_layout()
    plt.xlabel(expanded_names[x_field]);
    plt.ylabel("VBS speedup vs cuBLAS GEMM \n (higher is better)");
    
    #plt.ylim([0,1.5]);
    #plt.xlim([0,0.3]);
    plt.legend()
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")


dothis = True
if dothis:
    y_field = "A_blocks_density";
    y_field_function = lambda i: experiments[y_field][i];
    
    x_field = "A_total_nonzeros";
    x_field_function = mat_sparsity;
    
    plt.figure();
    A_cols = 4096; 
    A_rows = 4096;
    B_cols = 4096*4; 
    block_size = 256;
    accessibility = True
    plot_name = "scatter_performance_plot_" + str(block_size);
    if accessibility: plot_name += "_access";
    title = "\n" + "M,K,N = (" + str(A_rows) + "," + str(A_cols) + "," + str(B_cols) +")" + " \n Block size = " + str(block_size);

    fixed = {"B_cols": B_cols, "A_cols" : A_cols, "A_rows": A_rows, "A_block_size" : block_size};
    plot_best_scatter(x_field_function, y_field_function, "gemm_mean(ms)",fixed);
    plot_best_scatter(x_field_function, y_field_function, "cusparse_spmm_mean(ms)",fixed);
    plot_best_scatter(x_field_function, y_field_function, "VBSmm_mean(ms)",fixed);

    
    #drawing theoretical curve
    if accessibility: 
        x_s = np.linspace(0,0.1,100); 
        y_s = 1 - np.exp(-x_s*block_size)
        plt.plot(x_s,y_s)
        plt.fill_between(x_s, 1, y_s, alpha = 0.1, label = "accessible region")
    
    plt.tight_layout()
    plt.ylabel(expanded_names[y_field]);
    plt.xlabel("density of A");
    
    plt.ylim([0,1]);
    plt.xlim([0,0.1]);
    plt.legend()
    plt.title(title);
    if save:
        plt.savefig(saving_path + '/' + plot_name  + '.eps', format = 'eps', dpi=1000, bbox_inches = "tight")
        plt.savefig(saving_path + '/' + plot_name  + '.jpg', format = 'jpg', dpi=1000, bbox_inches = "tight")



