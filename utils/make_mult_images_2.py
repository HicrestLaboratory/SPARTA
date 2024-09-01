import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
import argparse

#----------------------ARGUMENTS
#______________________________________________________
parser = argparse.ArgumentParser(description="Analysis of multiplication times after reordering")
parser.add_argument("--routine", nargs="?", type=str, default="spmmcsr", help="the multiplication routine to report on, e.g. spmmcsr")
parser.add_argument("--root_dir", nargs="?", type=str, default="results/results_10-08-2024/", help="the directory where the csv of the experiments are stored")

args = parser.parse_args()

allowed_routines = ["spmmcsr", "spmvcsr", "spmmbsr", "spmvbsr"]
routine = args.routine
if routine not in allowed_routines:
    print(f"{routine} is not a valid routine {allowed_routines}")
    exit(1)
root_dir=args.root_dir

if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} does not exists.")
    exit(1)
#______________________________________________________

#____________PLOTTING PARAMS________________
color_dict = {
    "original": "blue",
    "clubs": "orange",
    "metis-edge-cut" : "red",
    "metis-volume" : "purple"
}

marker_dict = {
    "original": "o",
    "clubs": "x",
    "metis-edge-cut" : ".",
    "metis-volume" : "."
}




#----------------------FUNCTIONS
# _____________________________________________________
def count_best_method(best_dfs, methods):
    all_results_df = pd.DataFrame()
    common_matrices= find_common_matrices(best_dfs, methods)
    for method in methods:
        df = best_dfs[method].copy()
        df = set_allowed_matrices(df, common_matrices)
        df["method"] = method
        all_results_df = pd.concat([all_results_df,df])


    all_results_df = all_results_df.sort_values(by="method", ascending=False, key=lambda x: x == "original") #ensures ties are win by original
    overall_best_results = find_best(all_results_df)

    best_method_counts = {}
    for method in methods:
        best_method_counts[method] = (overall_best_results['method'] == method).sum()

    return best_method_counts

def calculate_speedups(best_dfs, methods):
    common_matrices = find_common_matrices(best_dfs,methods)
    speedups = {}
    for method in methods:
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df,common_matrices)
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        merged = pd.merge(original_df[['matrix', 'time']], method_df[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        merged['ratio'] = merged['time_method'] / merged['time_original']
        speedups[method] = 1./np.mean(merged['ratio'])        
    return speedups


def calculate_total_times(best_dfs, methods):
    sum_dict = {}
    common_matrices= find_common_matrices(best_dfs, methods)
    for method in methods:
        df = best_dfs[method].copy()
        clean_df = set_allowed_matrices(df,common_matrices)
        sum_dict[method] = sum(clean_df["time"])

    return sum_dict

def calculate_total_best_time(best_dfs, methods):
    common_matrices= find_common_matrices(best_dfs, methods)
    all_results_df = pd.DataFrame()
    for method in methods:
        df = best_dfs[method].copy()
        df = set_allowed_matrices(df,common_matrices)
        df["method"] =  method
        all_results_df = pd.concat([all_results_df,df])
    
    all_results_df = all_results_df.sort_values(by="method", ascending=False, key=lambda x: x == "original") #ensures ties are win by original
    overall_best_results = find_best(all_results_df)
    
    return sum(overall_best_results["time"]) 

# Helper function to read and concatenate files, ignoring empty ones
def read_and_concat(files):
    dfs = []
    for file in files:
        df = pd.read_csv(file, delim_whitespace=True, header=0)
        if not df.empty:
            dfs.append(df)
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    else:
        return pd.DataFrame()

def find_best(df):
    df = df.reset_index(drop=True)  # Ensure unique index
    return df.loc[df.groupby("matrix")["time"].idxmin()]

def set_allowed_matrices(df, allowed_matrices):
    return df[(df["matrix"]).isin(allowed_matrices)]

def find_common_matrices(dfs, methods):
    matrices_with_original = dfs["original"]['matrix'].unique()
    common_matrices=set(matrices_with_original)
    for method in methods:
        if method != "original":
            common_matrices &= set(dfs[method]['matrix'].unique())
    return common_matrices


def make_plot(values_dict, title="Plot", ylabel="Value", save_path = "test.png"):
    """    
    Parameters:
    values_dict (dict): A dictionary containing results for each method.
    """

    # Sort values by method names for consistent plotting
    methods_sorted = list(values_dict.keys())
    method_values = list(values_dict.values())

    # Plot the bar chart
    plt.figure(figsize=(10, 6))
    plt.bar(methods_sorted, method_values, color=[color_dict[method] for method in methods_sorted])
    
    # Add titles and labels
    plt.title(title)
    plt.xlabel("Reordering technique")
    plt.ylabel(ylabel)
    
    # Rotate the x-axis labels for better readability if necessary
    plt.xticks(rotation=45, ha="right")
    
    # Display the plot
    plt.tight_layout()
    plt.savefig(save_path)

def plot_speedup_distribution(best_dfs, methods, title = 'Distribution of Speedups Across Matrices', save_path= "test_dist.png"):
    """
    Plots the distribution (histogram) of speedups across matrices for each method.

    Parameters:
    - best_dfs: Dictionary of DataFrames for each method containing matrix times.
    - methods: List of methods to compare.
    - num_bins: Number of bins to use for the histogram (default is 10).
    """
    common_matrices = find_common_matrices(best_dfs, methods)
    plt.figure(figsize=(10, 6))
    
    for method in methods:
        if method == "original": continue
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df, common_matrices)
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        merged = pd.merge(original_df[['matrix', 'time']], method_df[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        merged['ratio'] = merged['time_method'] / merged['time_original']

        speedups = 1 / merged['ratio']
        speedups = [max(s, 1) for s in speedups]
        
        custom_bins = np.arange(1,3,0.05)
        plt.hist(speedups, bins=custom_bins, cumulative=True, color = color_dict[method], histtype='step', linewidth=2, label=f'{method}')
    plt.xlabel('Speedup')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)

def plot_speedup_vs_param(best_dfs, methods, param="density", title='Speedup vs matrix density', save_path="test_hist.png"):
    common_matrices = find_common_matrices(best_dfs, methods)
    plt.figure(figsize=(10, 6))
    
    for method in methods:
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df, common_matrices)
        
        # Include the 'param' column in the merge operation
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        
        # Merge including the 'param' column
        merged = pd.merge(original_df[['matrix', param, 'time']], 
                          method_df[['matrix', 'time']], 
                          on=['matrix'], 
                          suffixes=('_original', '_method'))
        
        merged = merged.sort_values(by=param)
        
        # Calculate the speedup ratio
        merged['ratio'] = merged['time_method'] / merged['time_original']
        speedups = 1 / merged['ratio']
        
        # Ensure speedups are capped at 1
        speedups[speedups < 1] = 1
        
        # Extract the values of the parameter to plot against
        param_values = merged[param]
        
        # Plot the data
        plt.plot(param_values, speedups, color = color_dict[method], marker='o', linestyle='',markersize= 3, label=f'{method}')
    
    plt.xscale('log')
    plt.ylim(0.9, 1.5)
    plt.xlabel(param.capitalize())
    plt.ylabel('Speedup')
    plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def plot_speedup_by_matrix(best_dfs, methods, title='Speedup for each matrix (sorted)', save_path="test_speedup_by_matrix.png"):
    """
    Plots the speedup for each matrix, sorted by speedups, for each method.

    Parameters:
    - best_dfs: Dictionary of DataFrames for each method containing matrix times.
    - methods: List of methods to compare.
    """
    common_matrices = find_common_matrices(best_dfs, methods)
    plt.figure(figsize=(10, 6))
    
    for method in methods:
        if method == "original" : continue
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df, common_matrices)
        
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        
        merged = pd.merge(original_df[['matrix', 'time']], 
                          method_df[['matrix', 'time']], 
                          on=['matrix'], 
                          suffixes=('_original', '_method'))
        
        # Calculate the speedup ratio
        merged['ratio'] = merged['time_method'] / merged['time_original']
        merged['speedup'] = 1 / merged['ratio']
        
        # Ensure speedups are capped at 1
        merged.loc[merged['speedup'] < 1, 'speedup'] = 1
        
        # Sort the DataFrame by speedup
        merged_sorted = merged.sort_values(by='speedup')
        
        # Plot the sorted speedup values
        plt.plot(range(len(merged_sorted)), merged_sorted['speedup'], color = color_dict[method], marker='o', linestyle='', markersize=3, label=f'{method}')
    
    plt.axhline(1, color = color_dict["original"])
    plt.ylim(0.9, 2)
    plt.yscale("log")
    plt.xlabel('Matrix (sorted by speedup)')
    plt.ylabel('Speedup')
    plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()





#-____________________ END OF FUNCTIONS_____________________________




methods=["original", "clubs", "metis-edge-cut", "metis-volume"]

#subdirectories
csv_dir = root_dir + "mult_csv/" + routine
output_plot_dir = root_dir + "mult_plots/" + routine
os.makedirs(output_plot_dir, exist_ok=True)

dfs = {method: pd.DataFrame() for method in methods}
files ={method: [os.path.join(csv_dir + "/" + method + "/", f) for f in os.listdir(csv_dir  + "/" + method + "/")] for method in methods}


#Fill all dataframes
for method in methods:
    dfs[method] = read_and_concat(files[method])
    dfs[method]["density"] = dfs[method]["nnz"]/(dfs[method]["rows"]*dfs[method]["cols"])

#metis only accepts square matrices
for method in ["metis-edge-cut", "metis-volume"]:
    dfs[method] = dfs[method][dfs[method]['rows'] == dfs[method]['cols']] 

#fix the tau for clubs 
#dfs["clubs"] = dfs["clubs"][dfs["clubs"]["tau"] == 1]


#clean data
cutoff = 3
for method in methods: 
    #remove unorderable matrices
    dfs[method] = dfs[method].loc[dfs[method]["nnz"] >= cutoff*dfs[method]["rows"]]

#find matrices for which original exist
matrices_with_original = dfs["original"]['matrix'].unique()
for method in methods:
    dfs[method] = dfs[method][(dfs[method]["matrix"]).isin(matrices_with_original)]
    allowed_matrices = len(dfs[method]["matrix"].unique())
    print(f"Method: {method}; Found data for {allowed_matrices} matrices")

#find matrices common among ALL methods
common_matrices=set(matrices_with_original)
for method in methods:
    if method != "original":
        common_matrices &= set(dfs[method]['matrix'].unique())
print(f"Found {len(common_matrices)} common matrices")


#find matrices that exist only for a method
method_matrices = {method : set(dfs[method]["matrix"].unique()) - common_matrices for method in methods}
best_dfs = {}
for method in methods:
    best_dfs[method] = find_best(dfs[method])



comparisons = [methods, ["original", "clubs", "metis-edge-cut"], ["original", "clubs"], ["original", "metis-edge-cut", "metis-volume"]]

for exp_methods in comparisons:
    exp_methods_string = "_".join(exp_methods)
    print("*****************")
    print("Comparing: ", exp_methods)

    n_matrices = len(find_common_matrices(best_dfs, exp_methods))
    print(f"Found {n_matrices} common matrices")

    counts = count_best_method(best_dfs, exp_methods)
    print("BEST COUNT ON COMMON MATRICES: ", counts)
    make_plot(values_dict=counts, 
              title=f"Best reordering method for {routine} on {n_matrices} matrices for {exp_methods}",
              ylabel="# of times X is the best method",
              save_path= f"{output_plot_dir}/{routine}_best_count_plot_{exp_methods_string}")


    speedups = calculate_speedups(best_dfs, exp_methods)
    print("SPEEDUPS GEOMETRIC MEAN ON COMMON MATRICES: ", speedups)
    make_plot(values_dict=speedups, 
              title=f"Geometric mean of {routine} speedup on {n_matrices} matrices for {exp_methods}",
              ylabel="Speed-up against non-reordered matrix",
              save_path= f"{output_plot_dir}/{routine}_speedup_geomean_{exp_methods_string}")

    total_times = calculate_total_times(best_dfs, exp_methods)
    print("TOTAL TIMES ON COMMON MATRICES:", total_times)
    make_plot(values_dict=total_times, 
              title=f"Total time to run {routine} on {n_matrices} matrices for {exp_methods}",
              ylabel=f"Total {routine} time",
              save_path= f"{output_plot_dir}/{routine}_total_times_{exp_methods_string}")

    print("TOTAL TIME USING BEST: ", calculate_total_best_time(best_dfs, exp_methods))

    plot_speedup_distribution(best_dfs, 
                              exp_methods, 
                              title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",

                              save_path= f"{output_plot_dir}/{routine}_cumulative_hist_{exp_methods_string}.png"
                              )

    param = "rows"
    plot_speedup_vs_param(  best_dfs, 
                            exp_methods, 
                            param=param,
                            title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_speedup_vs_{param}_{exp_methods_string}.png"
                              )

    param = "density"
    plot_speedup_vs_param(  best_dfs, 
                            exp_methods, 
                            param=param,
                            title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_speedup_vs_{param}_{exp_methods_string}.png"
                              )
    
    plot_speedup_by_matrix( best_dfs, 
                            exp_methods, 
                            title=f"Speedup by matrix for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_speedup_by_matrix_{exp_methods_string}.png"
                              )




#TRY: CALCULATE TOTAL INCLUDING RECTANGULAR FOR METIS










