import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean

routine="spmmcsr"
methods=["original", "clubs", "saad", "metis-edge-cut", "metis-volume"]

# Define the directory containing the files
directory = 'results/results_2024/mult/' + routine
output_plot_dir = "results/results_2024/plots/mult/" + routine
os.makedirs(output_plot_dir, exist_ok=True)


columns={}
columns["clubs"] = ["routine", "matrix", "algo", "mask", "centroid", "tau", "rows", "cols", "nnz", "time"]
columns["original"] = ["routine", "matrix", "algo", "rows", "cols", "nnz", "time"]
columns["saad"] = ["routine", "matrix", "algo", "tau", "rows", "cols", "nnz", "time"]
columns["metis-volume"] = ["routine", "matrix", "algo", "objective", "parts", "rows", "cols", "nnz", "time"]
columns["metis-edge-cut"] = ["routine", "matrix", "algo", "objective", "parts", "rows", "cols", "nnz", "time"]

dfs = {method: pd.DataFrame() for method in methods}
files ={method: [os.path.join(directory + "/" + method + "/", f) for f in os.listdir(directory  + "/" + method + "/")] for method in methods}

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

#Fill all dataframes
for method in methods:
    dfs[method] = read_and_concat(files[method])
    print(dfs[method])

# Function to get the best result for each matrix
def get_best_results(df, group_cols, value_col='time'):
    return df.loc[df.groupby(group_cols)[value_col].idxmin()]

# Calculate the geometric mean of time ratios
def calculate_geometric_mean_ratios(dfs, methods):
    ratios = {method: [] for method in methods if method != 'original'}
    original_best = get_best_results(dfs['original'], ['matrix'])

    for method in methods:
        if dfs[method].empty:
            continue
        if method == 'original':
            continue
        best_results = get_best_results(dfs[method],'matrix')
        merged = pd.merge(original_best[['matrix', 'time']], best_results[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        merged['ratio'] = merged['time_method'] / merged['time_original']
        #merged = merged[merged["ratio"] < 1] #only counts effective reorderings
        #merged["ratio"][merged["ratio"] > 1] = 1 #take the original when reordering is bad
        ratios[method] = gmean(merged['ratio'])

    return ratios

# Function to count how many matrices each method is the best
def count_best_methods(dfs, methods):
    best_method_counts = {method: 0 for method in methods}
    all_best_results = pd.DataFrame()

    for method in methods:
        if not dfs[method].empty:
            best_results = get_best_results(dfs[method], ['matrix'])
            best_results['method'] = method
            all_best_results = pd.concat([all_best_results, best_results])

    overall_best_results = get_best_results(all_best_results, ['matrix'], 'time')
    
    for method in methods:
        best_method_counts[method] = (overall_best_results['method'] == method).sum()

    return best_method_counts

# Function to count unique matrices for each method
def count_unique_matrices(dfs, methods):
    unique_matrix_counts = {}
    for method in methods:
        if not dfs[method].empty:
            unique_matrix_counts[method] = dfs[method]['matrix'].nunique()
        else:
            unique_matrix_counts[method] = 0
    return unique_matrix_counts

# Calculate geometric mean ratios
geometric_mean_ratios = calculate_geometric_mean_ratios(dfs, methods)
print(geometric_mean_ratios)

# Count best methods
best_method_counts = count_best_methods(dfs, methods)
print(best_method_counts)

# Count unique matrices for each method
unique_matrix_counts = count_unique_matrices(dfs, methods)
print("Unique matrix counts:", unique_matrix_counts)






