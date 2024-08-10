import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean

routine="spmvcsr"
methods=["original", "clubs", "saad", "metis-edge-cut", "metis-volume", "denseAMP"]

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

#metis only accepts square matrices
for method in ["metis-edge-cut", "metis-volume"]:
    dfs[method] = dfs[method][dfs[method]['rows'] == dfs[method]['cols']] 

#remove matrices for which original does not exist
matrices_with_original = dfs["original"]['matrix'].unique()
for method in methods:
    dfs[method] = dfs[method][(dfs[method]["matrix"]).isin(matrices_with_original)]
    allowed_matrices = len(dfs[method]["matrix"].unique())
    print(f"Method: {method}; Found data for {allowed_matrices} matrices")


def find_best(df):
    return df.loc[df.groupby("matrix")["time"].idxmin()]
    
    
#TRY: CALCULATE TOTAL INCLUDING RECTANGULAR FOR METIS
    