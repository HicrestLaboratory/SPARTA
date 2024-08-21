import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean

def count_best_method(best_dfs, matrices):
    all_results_df = pd.DataFrame()
    for method in methods:
        df = best_dfs[method].copy()
        df = set_allowed_matrices(df, matrices)
        df["method"] = method
        all_results_df = pd.concat([all_results_df,df])


    all_results_df = all_results_df.sort_values(by="method", ascending=False, key=lambda x: x == "original") #ensures ties are win by original
    overall_best_results = find_best(all_results_df)

    best_method_counts = {}
    for method in methods:
        best_method_counts[method] = (overall_best_results['method'] == method).sum()

    return best_method_counts

def calculate_total_times(best_dfs, matrices):
    sum_dict = {}
    for method in methods:
        df = best_dfs[method].copy()
        clean_df = set_allowed_matrices(df,matrices)
        sum_dict[method] = sum(clean_df["time"])

    return sum_dict

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

#-____________________ END OF FUNCTIONS_____________________________




routine="spmvcsr"
methods=["original", "clubs", "metis-edge-cut", "metis-volume"]

# Define the directory containing the files
root_dir="results/results_10-08-2024/"
if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} exists.")
    exit(1)

csv_dir = root_dir + "mult_csv/" + routine
output_plot_dir = root_dir + "mult_plots/" + routine
os.makedirs(output_plot_dir, exist_ok=True)

dfs = {method: pd.DataFrame() for method in methods}
files ={method: [os.path.join(csv_dir + "/" + method + "/", f) for f in os.listdir(csv_dir  + "/" + method + "/")] for method in methods}


#Fill all dataframes
for method in methods:
    dfs[method] = read_and_concat(files[method])

#metis only accepts square matrices
for method in ["metis-edge-cut", "metis-volume"]:
    dfs[method] = dfs[method][dfs[method]['rows'] == dfs[method]['cols']] 

#clean data
cutoff = 2
for method in methods:
    
    #remove diagonal matrices
    dfs[method] = dfs[method].loc[dfs[method]["nnz"] >= 2*dfs[method]["rows"]]

    #remove long
    #dfs[method] = dfs[method][dfs[method]["time"] < 0.05]

#remove matrices for which original does not exist
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

method_matrices = {method : set(dfs[method]["matrix"].unique()) - common_matrices for method in methods}

best_dfs = {}
for method in methods:
    best_dfs[method] = find_best(dfs[method])


        


print("ALL BEST: ", count_best_method(best_dfs, matrices_with_original))

print("COMMON BEST: ", count_best_method(best_dfs, common_matrices))
print("TOTAL TIMES COMMON:", calculate_total_times(best_dfs, common_matrices))

print("CLUBS MATRICES: ", count_best_method(best_dfs, method_matrices["clubs"]))
print("TOTAL TIMES CLUBS:", calculate_total_times(best_dfs, method_matrices["clubs"]))


#TRY: CALCULATE TOTAL INCLUDING RECTANGULAR FOR METIS










