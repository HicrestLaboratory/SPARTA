import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean

methods=["original", "clubs", "saad", "metis-edge-cut", "metis-volume"]

# Define the directory containing the files
directory = 'results/results_2024/reorder_csv_25_07'
output_plot_dir = 'results/results_2024/plots'
os.makedirs(output_plot_dir, exist_ok=True)

dfs = {method: pd.DataFrame() for method in methods}
files = {}
for method in ["original", "clubs", "saad"]:
    files[method] = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith(method)]
for method in ["metis-edge-cut", "metis-volume"]:
    files[method] = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith("metis")]


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
    print(method, dfs[method])

#remove mismatched metis objectives
for metis_obj in ["edge-cut","volume"]:
    key = "metis-" + metis_obj
    df = dfs[key]
    dfs[key] = df[df['metis_obj'] != metis_obj]

# Function to get the best result for each matrix
def get_best_results(df, group_cols, value_col='VBR_nzblocks_count'):
    return df.loc[df.groupby(group_cols)[value_col].idxmin()]

# Calculate the geometric mean of time ratios
def calculate_geometric_mean_ratios(dfs, methods, value_col="VBR_nzblocks_count"):
    ratios = {method: [] for method in methods if method != 'original'}
    original_best = get_best_results(dfs['original'], ['matrix_name'])

    for method in methods:
        if dfs[method].empty:
            continue
        if method == 'original':
            continue
        best_results = get_best_results(dfs[method],'matrix_name', value_col=value_col)
        merged = pd.merge(original_best[['matrix_name', value_col]], best_results[['matrix_name', value_col]], on='matrix_name', suffixes=('_original', '_method'))
        merged['ratio'] = merged[value_col + '_method'] / merged[value_col + '_original']
        #merged = merged[merged["ratio"] < 1] #only counts effective reorderings
        #merged["ratio"][merged["ratio"] > 1] = 1 #take the original when reordering is bad
        ratios[method] = gmean(merged['ratio'])

    return ratios

# Function to count how many matrices each method is the best
def count_best_methods(dfs, methods, value_col="VBR_nzblocks_count"):
    best_method_counts = {method: 0 for method in methods}
    all_best_results = pd.DataFrame()

    for method in methods:
        if not dfs[method].empty:
            best_results = get_best_results(dfs[method], ['matrix_name'])
            best_results['method'] = method
            all_best_results = pd.concat([all_best_results, best_results])

    overall_best_results = get_best_results(all_best_results, ['matrix_name'], value_col)
    
    for method in methods:
        best_method_counts[method] = (overall_best_results['method'] == method).sum()

    return best_method_counts

# Function to count unique matrices for each method
def count_unique_matrices(dfs, methods):
    unique_matrix_counts = {}
    for method in methods:
        if not dfs[method].empty:
            unique_matrix_counts[method] = dfs[method]['matrix_name'].nunique()
        else:
            unique_matrix_counts[method] = 0
    return unique_matrix_counts

# Calculate geometric mean ratios
print("GEOMETRIC MEANS:")
for value_col in ["VBR_nzblocks_count","VBR_longest_row"]:
    geometric_mean_ratios = calculate_geometric_mean_ratios(dfs, methods, value_col=value_col)
    print(value_col + ":", geometric_mean_ratios)
   

print("BEST METHOD MATRIX COUNT")
# Count best methods
for value_col in ["VBR_nzblocks_count","VBR_longest_row"]:
    best_method_counts = count_best_methods(dfs, methods, value_col)
    print(value_col + ":", best_method_counts)

# Count unique matrices for each method
unique_matrix_counts = count_unique_matrices(dfs, methods)
print("Unique matrix counts:", unique_matrix_counts)






