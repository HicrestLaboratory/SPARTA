import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean

routine="spmmcsr"
methods=["original", "clubs", "metis-edge-cut", "metis-volume", "saad"]
#methods=["original", "clubs", "saad", "metis-edge-cut", "metis-volume", "denseAMP"]

# Define the directory containing the files
root_dir="results/results_10-08-2024/"
if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} exists.")
    exit(1)

csv_dir = root_dir + "mult_csv/" + routine
output_plot_dir = root_dir + "mult_plots/" + routine
os.makedirs(output_plot_dir, exist_ok=True)


columns={}
columns["clubs"] = ["routine", "matrix", "algo", "mask", "centroid", "tau", "rows", "cols", "nnz", "time"]
columns["original"] = ["routine", "matrix", "algo", "rows", "cols", "nnz", "time"]
columns["saad"] = ["routine", "matrix", "algo", "tau", "rows", "cols", "nnz", "time"]
columns["metis-volume"] = ["routine", "matrix", "algo", "objective", "parts", "rows", "cols", "nnz", "time"]
columns["metis-edge-cut"] = ["routine", "matrix", "algo", "objective", "parts", "rows", "cols", "nnz", "time"]

dfs = {method: pd.DataFrame() for method in methods}
files ={method: [os.path.join(csv_dir + "/" + method + "/", f) for f in os.listdir(csv_dir  + "/" + method + "/")] for method in methods}

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
    #print(dfs[method])


#metis only accepts square matrices
for method in ["metis-edge-cut", "metis-volume"]:
    if method in dfs.keys():
        dfs[method] = dfs[method][dfs[method]['rows'] == dfs[method]['cols']] 

#cut diagonal matrices
for method in methods:
    dfs[method] = dfs[method].loc[dfs[method]["nnz"] >= 16*dfs[method]["rows"]]

# Function to get the best result for each matrix
def get_best_results(df, group_cols, value_col='time'):
    return df.loc[df.groupby(group_cols)[value_col].idxmin()]

# Calculate the geometric mean of time ratios
def calculate_geometric_mean_ratios(dfs, methods, common: bool = False):
    ratios = {method: [] for method in methods if method != 'original'}
    original_best = get_best_results(dfs['original'], ['matrix'])

    common_matrices = set(dfs[methods[0]]['matrix'].unique())
    for method in methods[1:]:
        common_matrices &= set(dfs[method]['matrix'].unique())

    for method in methods:
        if dfs[method].empty:
            continue
        best_results = get_best_results(dfs[method],'matrix')
        if common: best_results = best_results[best_results['matrix'].isin(common_matrices)]
        merged = pd.merge(original_best[['matrix', 'time']], best_results[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        #print("MERGED DATASET", method, merged)
        merged['ratio'] = merged['time_method'] / merged['time_original']
        #merged = merged[merged["ratio"] <= 1] #only counts effective reorderings
        #merged.loc[merged["ratio"] > 1, "ratio"] = 1
        ratios[method] = np.mean(merged['ratio'])

    return ratios

# Calculate the total time
def calculate_total_time_ratio(dfs, methods, common : bool = False):
    ratios = {method: [] for method in methods if method != 'original'}
    original_best = get_best_results(dfs['original'], ['matrix'])

    common_matrices = set(dfs[methods[0]]['matrix'].unique())
    for method in methods[1:]:
        common_matrices &= set(dfs[method]['matrix'].unique())

    for method in methods:
        if dfs[method].empty:
            continue
        best_results = get_best_results(dfs[method],'matrix')
        if common: best_results = best_results[best_results['matrix'].isin(common_matrices)]
        merged = pd.merge(original_best[['matrix', 'time']], best_results[['matrix', 'time']], on='matrix', suffixes=('_original', '_method'))
        #print("MERGED DATASET", method, merged)
        merged['ratio'] = merged['time_method'] / merged['time_original']
        #merged = merged[merged["ratio"] <= 1] #only counts effective reorderings
        #merged.loc[merged["ratio"] > 1, "ratio"] = 1
        ratios[method] = sum(merged['time_method'])/sum(merged["time_original"])

    return ratios


# Function to count how many matrices each method is the best
def count_best_methods(dfs, methods, common : bool = False):
    best_method_counts = {method: 0 for method in methods}
    all_best_results = pd.DataFrame()

    common_matrices = set(dfs[methods[0]]['matrix'].unique())
    for method in methods[1:]:
        common_matrices &= set(dfs[method]['matrix'].unique())

    for method in methods:
        if not dfs[method].empty:
            best_results = get_best_results(dfs[method], ['matrix'])
            if common: best_results = best_results[best_results['matrix'].isin(common_matrices)]
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


# Count unique matrices for each method
unique_matrix_counts = count_unique_matrices(dfs, methods)
print("Unique matrix counts:", unique_matrix_counts)


matrices_with_original = dfs["original"]['matrix'].unique()
print("ALLOWED MATRICES (original time exists): ", len(matrices_with_original))



for common in [True, False]:
    if common: print("************Calculating only on common matrices")
    else: print("************Calculating on all matrices")

    # Calculate geometric mean ratios
    geometric_mean_ratios = calculate_geometric_mean_ratios(dfs, methods, common)
    print("geometric_mean_ratios", geometric_mean_ratios)

    total_time_ratios = calculate_total_time_ratio(dfs, methods, common)
    print("total_time_ratios", total_time_ratios)


    #Count best methods
    best_method_counts = count_best_methods(dfs, methods, common)
    print("best_method_counts", best_method_counts)






