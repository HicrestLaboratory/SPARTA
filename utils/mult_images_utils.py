import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
import argparse
import seaborn as sns
from matplotlib.ticker import FuncFormatter

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
        merged['ratio'] = merged['time_original'] / merged['time_method']
        speedups[method] = gmean(merged['ratio'].values)  
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

def find_best(df, params = []):
    df = df.reset_index(drop=True)  # Ensure unique index
    params = ["matrix",] + params
    print(params)
    return df.loc[df.groupby(params)["time"].idxmin()]

def set_allowed_matrices(df, allowed_matrices):
    return df[(df["matrix"]).isin(allowed_matrices)]

def find_common_matrices(dfs, methods):
    matrices_with_original = dfs["original"]['matrix'].unique()
    common_matrices=set(matrices_with_original)
    for method in methods:
        if method != "original":
            common_matrices &= set(dfs[method]['matrix'].unique())
    return common_matrices


def make_plot(values_dict, title="Plot", ylabel="Value", percent = False, save_path = "test.png"):
    """    
    Parameters:
    values_dict (dict): A dictionary containing results for each method.
    """

    # Sort values by method names for consistent plotting
    methods_sorted = list(values_dict.keys())
    method_values = list(values_dict.values())

    # Plot the bar chart
    plt.figure(figsize=(10, 6))
    
    # Add titles and labels
    #plt.title(title)
    plt.xlabel("Reordering technique")
    plt.ylabel(ylabel)
    
    if percent:
        method_values = [m*100 - 100 for m in method_values]
        plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    # Rotate the x-axis labels for better readability if necessary

    
    plt.bar(methods_sorted, method_values, color=[color_dict[method] for method in methods_sorted])
    #plt.xticks(rotation=45, ha="right")
    
    # Display the plot
    plt.tight_layout()
    plt.savefig(save_path)

def make_speedup_barplot(best_dfs, methods, title="Plot", ylabel="Value", save_path = "test.png"):
    """    
    Parameters:
    values_dict (dict): A dictionary containing results for each method.
    """

    common_matrices = find_common_matrices(best_dfs,methods)
    speedups = {}
    low_errors = []
    up_errors = []
    for method in methods:
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df,common_matrices)
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        merged = pd.merge(original_df[['matrix', 'time']], method_df[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        merged['ratio'] = merged['time_original'] / merged['time_method']
        #merged['ratio'] = merged['ratio']*100 - 100
        percentiles = np.percentile(merged['ratio'].values, [25, 50, 75])
        #speedups[method] = gmean(merged['ratio']) 
        speedups[method] = percentiles[1] 
        print(method, percentiles, speedups[method])
        low_errors.append(speedups[method] - percentiles[0])
        up_errors.append(percentiles[2] - speedups[method])
       
    # Plot the bar chart
    plt.figure(figsize=(10, 6))

    # Add titles and labels
    #plt.title(title)
    plt.xlabel("Reordering technique")
    plt.ylabel(ylabel)
    
    plt.bar(methods, speedups.values(), yerr=[low_errors, up_errors], color=[color_dict[method] for method in methods])
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y*100 - 100:.0f}%'))
    #plt.xticks(rotation=45, ha="right")
    max_y = (max(up_errors) + max(speedups.values()))*1.01
    plt.ylim(1,max_y)
    # Display the plot
    plt.tight_layout()
    plt.savefig(save_path)

def make_speedup_violin(best_dfs, methods, title="Plot", ylabel="Value", save_path = "test.png"):
    """ 
    Parameters:
    best_dfs (dict): A dictionary containing results for each method.
    methods (list): A list of methods to compare.
    title (str): Title of the plot.
    ylabel (str): Y-axis label.
    save_path (str): Path to save the plot.
    """
    common_matrices = find_common_matrices(best_dfs, methods)
    speedup_data = []
    method_labels = []
    
    for method in methods:
        # Preprocess the data
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df, common_matrices)
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        merged = pd.merge(original_df[['matrix', 'time']], method_df[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        
        # Calculate the speedup (original time / method time)
        merged['ratio'] = merged['time_original'] / merged['time_method']
        merged.loc[merged["ratio"] < 1] = 1

        # Collect data for the violin plot
        speedup_data.extend(merged['ratio'].values)
        method_labels.extend([method] * len(merged))

    # Create a DataFrame for seaborn
    plot_df = pd.DataFrame({'Method': method_labels, 'Speedup': speedup_data})

    # Plot the violin plot
    plt.figure(figsize=(10, 6))
    ax = sns.violinplot(x='Method', y='Speedup', data=plot_df, palette=[color_dict[method] for method in methods], )
    plt.setp(ax.collections, alpha=.3)


    # Add titles and labels
    #plt.title(title)
    plt.xlabel("Reordering Technique")
    plt.ylabel(ylabel)

    # Format the y-axis to show percentages (if needed)
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y*100 - 100:.0f}%'))

    # Adjust layout and save the plot
    #plt.xticks(rotation=45, ha="right")
    plt.ylim(0.8,1.5)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_speedup_distribution(best_dfs, methods, matrices = None, title = 'Distribution of Speedups Across Matrices', save_path= "test_dist.png"):
    """
    Plots the distribution (histogram) of speedups across matrices for each method.

    Parameters:
    - best_dfs: Dictionary of DataFrames for each method containing matrix times.
    - methods: List of methods to compare.
    - num_bins: Number of bins to use for the histogram (default is 10).
    """
    common_matrices = find_common_matrices(best_dfs, methods)
    if matrices != None: 
        common_matrices &= matrices

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
        speedups_percent = [s*100 - 100 for s in speedups]
        
        custom_bins = np.arange(0,200,5)
        plt.hist(speedups_percent, bins=custom_bins, cumulative=False, color = color_dict[method], histtype='step', linewidth=2, label=f'{method}')
    
    plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    plt.xlabel('Speedup')
    plt.ylabel('Frequency')
    #plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)

def plot_speedup_vs_matrix_param(best_dfs, methods, matrices = None, param="density", title='Speedup vs matrix density', save_path="test_hist.png"):

    common_matrices = find_common_matrices(best_dfs, methods)
    if matrices != None: 
        common_matrices &= matrices

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
    #plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()


def plot_speedup_by_matrix(best_dfs, order_by, methods, matrices = None, title='Speedup for each matrix (sorted)', save_path="test_speedup_by_matrix.png"):
    """
    Plots the speedup for each matrix, sorted by speedups, for each method.

    Parameters:
    - best_dfs: Dictionary of DataFrames for each method containing matrix times.
    - methods: List of methods to compare.
    """
    common_matrices = find_common_matrices(best_dfs, methods)
    if matrices != None: 
        common_matrices &= matrices
    
    # First, calculate the speedup for the 'clubs' method and store the ordering
    if not order_by in methods: 
        print("WARNING: plot_speedup_by_matrix: ORDER_BY must appear also in METHODS")
        return
    ordered_df = best_dfs[order_by].copy()
    ordered_df = set_allowed_matrices(ordered_df, common_matrices)
    
    original_df = best_dfs["original"].copy()
    original_df = set_allowed_matrices(original_df, common_matrices)
    
    merged_ordered_by = pd.merge(original_df[['matrix', 'time']], 
                            ordered_df[['matrix', 'time']], 
                            on=['matrix'], 
                            suffixes=('_original', '_method'))
    
    # Calculate the speedup ratio for 'clubs'
    merged_ordered_by['ratio'] = merged_ordered_by['time_method'] / merged_ordered_by['time_original']
    merged_ordered_by['speedup'] = 1 / merged_ordered_by['ratio']
    
    # Ensure speedups are capped at 1
    #merged_ordered_by.loc[merged_ordered_by['speedup'] < 1, 'speedup'] = 1
    
    # Sort the DataFrame by speedup for 'clubs'
    merged_clubs_sorted = merged_ordered_by.sort_values(by='speedup')
    ordered_matrices = merged_clubs_sorted['matrix'].values  # Store the ordered list of matrices
    
    plt.figure(figsize=(10, 6))
    # Now, plot the speedups for all methods using the same ordering
    for method in methods:
        if method == "original": 
            continue
        
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        
        merged = pd.merge(original_df[['matrix', 'time']], 
                          method_df[['matrix', 'time']], 
                          on=['matrix'], 
                          suffixes=('_original', '_method'))
        
        # Calculate the speedup ratio
        merged['ratio'] = merged['time_method'] / merged['time_original']
        merged['speedup'] = 1 / merged['ratio']
        merged['speedup_percent'] = merged['speedup'] * 100 - 100

        # Ensure speedups are capped at 1
        #merged.loc[merged['speedup'] < 1, 'speedup'] = 1
        
        # Reorder the DataFrame to match the 'clubs' order
        merged['matrix'] = pd.Categorical(merged['matrix'], categories=ordered_matrices, ordered=True)
        merged_sorted = merged.sort_values(by='matrix')
        
        # Plot the sorted speedup values
        plt.plot(range(len(merged_sorted)), merged_sorted['speedup_percent'], color=color_dict[method], marker='o', linestyle='', markersize=3, label=f'{method}')
    
    plt.axhline(0, color=color_dict["original"])
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    plt.ylim(-10, 50)
    plt.xlabel(f'Matrix (sorted by speedup on {order_by})')
    plt.ylabel('Speedup')
    #plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_club_performance_by_param(dfs, method, param = "mask", title='Performance Variation under Mask', save_name=""):
    """
    Plots the performance variation for a specific method under different mask values,
    with one curve for each tau value, keeping centroid fixed.
    """
    # Filter the DataFrame for the selected method and fixed centroid value
    
    method_df = dfs[method].copy()
    original_df = dfs["original"]
    method_df = find_best(method_df, [param,])
    
    plt.figure(figsize=(10, 6))

    merged = pd.merge(original_df[['matrix', 'time']], 
                        method_df[['matrix', param,'time']], 
                        on=['matrix'], 
                        suffixes=('_original', '_method'))
    
    merged = merged.sort_values(by=param)

    merged['speedup'] = np.maximum(1, merged['time_original'] / merged['time_method'])

    # Ensure speedups are capped at 1

    plot_data = merged[[param,"speedup"]]

    #adds best plot
    best_merged = merged.loc[merged.groupby(by="matrix")["speedup"].idxmax()]
    speeds = best_merged["speedup"].values
    best_data = pd.DataFrame({param: ["Best"] * len(best_merged), "speedup": best_merged["speedup"].values})
    #---------------

    plot_data = pd.concat([plot_data, best_data], ignore_index=True)
    plot_data['speedup_percent'] = plot_data['speedup'] * 100 - 100

    sns.violinplot(x=param, y='speedup_percent', data=plot_data, inner='quartile', color=color_dict["clubs"])
# Plot the data with asymmetric error bars
    #plt.errorbar(param_values, speedups, yerr=[lower_errors, upper_errors], color=color_dict["clubs"], marker='o', linestyle='--', markersize=3, capsize=5)
    
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    plt.ylim(-10,50)
    plt.xlabel(param)
    plt.ylabel('Geometric mean of the speedup (%)')
    #plt.title(title)
    plt.tight_layout()
    plt.savefig(save_name)
    plt.close()


    return 
    #lineplot
    speedups = []
    lower_errors = []
    upper_errors = []

    param_values = merged[param].values
    for param_value in param_values:
        speeds = merged[merged[param] == param_value]["speedup"].values
        speedup = gmean(speeds)
        speedups.append(speedup)
        
        # Calculate percentiles for asymmetric error bars
        lower_error = speedup - np.percentile(speeds, 25)
        upper_error = np.percentile(speeds, 75) - speedup
        # Ensure errors are non-negative
        lower_errors.append(max(0, lower_error))
        upper_errors.append(max(0, upper_error))

#-____________________ END OF FUNCTIONS_____________________________

