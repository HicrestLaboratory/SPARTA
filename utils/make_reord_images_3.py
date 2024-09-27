import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
import argparse
from images_utils import *

#----------------------ARGUMENTS
#______________________________________________________
parser = argparse.ArgumentParser(description="Analysis of multiplication times after reordering")
parser.add_argument("--root_dir", nargs="?", type=str, default="results/results_10-08-2024/", help="the directory where the csv of the experiments are stored")
parser.add_argument("--bsize", nargs="?", type=str, default="64", help="size of the dense matrix")

args = parser.parse_args()

root_dir=args.root_dir
bsize=args.bsize

if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} does not exists.")
    exit(1)
#______________________________________________________


methods=["original", "clubs", "metis-edge-cut", "metis-volume"]

#subdirectories
csv_reordering_dir = root_dir + "reorder_csv"
csv_multiplication_dir = root_dir + "mult_csv"
output_plot_dir = root_dir + "reorder_plots/"
os.makedirs(output_plot_dir, exist_ok=True)

dfs_reordering = {method: pd.DataFrame() for method in methods}
dfs_multiplication = {method: pd.DataFrame() for method in methods}

#Fill all dataframes
for method in methods:
    reordering_file = f"{csv_reordering_dir}/{method}_scramble0_bsize{bsize}.txt"
    dfs_reordering[method] = pd.read_csv(reordering_file, delim_whitespace=True, header=0)
    dfs_reordering[method].rename(columns={'matrix_name': 'matrix'}, inplace=True)
    
    multiplication_file = f"{csv_multiplication_dir}/spmmbsr/{method}/" + os.listdir(csv_multiplication_dir  + "/spmmbsr/" + method + "/")[0]
    dfs_multiplication[method] = pd.read_csv(multiplication_file, delim_whitespace=True, header=0)


#add the speedups
for method in methods:


    common_cols = set(dfs_reordering[method].columns) & set(dfs_multiplication[method].columns)
    dfs_multiplication[method] = dfs_multiplication[method].merge(
        dfs_multiplication["original"][["matrix",'time']],
        on="matrix",
        how="left",
        suffixes=('', '_original')
    )

    print("dfs_reordering[method] merge keys:")
    print(dfs_reordering[method][list(common_cols)].drop_duplicates())

    print("dfs_multiplication[method] merge keys:")
    print(dfs_multiplication[method][list(common_cols)].drop_duplicates())

    dfs_multiplication[method]['speedup'] = dfs_multiplication[method]['time_original'] / dfs_multiplication[method]['time']
    print(common_cols)
    dfs_reordering[method] = dfs_reordering[method].merge(
        dfs_multiplication[method][list(common_cols) + ["speedup"]],
        on=list(common_cols),
        how="left"
    )
    print(dfs_reordering[method].columns)
    print(method, dfs_reordering[method]["speedup"].unique())


exit()



#metis only accepts square matrices
for method in ["metis-edge-cut", "metis-volume"]:
    dfs_reordering[method] = dfs_reordering[method][dfs_reordering[method]['rows'] == dfs_reordering[method]['cols']] 

#fix the tau for clubs 
#dfs["clubs"] = dfs["clubs"][dfs["clubs"]["tau"] == 1]

#clean data
degree_cutoff = 2
nnz_cutoff = 0
for method in methods: 
    #remove unorderable matrices
    dfs_reordering[method] = dfs_reordering[method].loc[dfs_reordering[method]["nnz"] >= degree_cutoff*dfs_reordering[method]["rows"]]
    dfs_reordering[method] = dfs_reordering[method].loc[dfs_reordering[method]["nnz"] >= nnz_cutoff]


all_matrices = set()
for method in methods:
    all_matrices = all_matrices.union(set(dfs_reordering[method]["matrix"].unique()))
print(f"TOTAL MATRIX COUNT: {len(all_matrices)}")


#find matrices for which original exist
matrices_with_original = dfs_reordering["original"]['matrix'].unique()
for method in methods:
    #dfs[method] = dfs[method][(dfs[method]["matrix"]).isin(matrices_with_original)]
    allowed_matrices = len(dfs_reordering[method]["matrix"].unique())
    print(f"Method: {method}; Found data for {allowed_matrices} matrices")

#find matrices common among ALL methods
common_matrices=set(matrices_with_original)
for method in methods:
    if method != "original":
        common_matrices &= set(dfs_reordering[method]['matrix'].unique())
print(f"Found {len(common_matrices)} common matrices")

#find rectangular and square matrices

rectangular_matrices = set()
square_matrices = set()
for method in methods: 
    rectangular_matrices |= set(dfs_reordering[method][dfs_reordering[method]['rows'] != dfs_reordering[method]['cols']]["matrix"].unique())
    square_matrices |= set(dfs_reordering[method][dfs_reordering[method]['rows'] == dfs_reordering[method]['cols']]["matrix"].unique())
print(f"RECTANGULAR: {len(rectangular_matrices)}, SQUARE: {len(square_matrices)}")


failed_matrices = {}
for method in methods:
    valid_matrices = set(all_matrices)
    working_matrices = set(dfs_reordering[method]["matrix"].unique()) 
    if "metis" in method:
        valid_matrices -= rectangular_matrices
        working_matrices -= rectangular_matrices
    failed_matrices[method] = valid_matrices - working_matrices
    print(f"METHOD: {method}, INVALID: {len(all_matrices - valid_matrices)}, FAILED: {len(failed_matrices[method])}, SUCCESS = {len(working_matrices)}")    
    

#find matrices that exist only for a method
method_matrices = {method : set(dfs_reordering[method]["matrix"].unique()) - common_matrices for method in methods}
best_dfs = {}
for method in methods:
    best_dfs[method] = find_best(dfs_reordering[method], best_parameter= "VBR_nzblocks_count")

comparisons = [methods, ["original", "clubs", "metis-edge-cut"], ["original", "clubs"], ["original", "metis-edge-cut", "metis-volume"]]


single_methods = [[m] for m in methods]
for exp_methods in single_methods:
    n_matrices = len(find_common_matrices(best_dfs, exp_methods))
    print(f"Found {n_matrices} common matrices")


for exp_methods in comparisons:
    exp_methods_string = "_".join(exp_methods)
    print("*****************")
    print("Comparing: ", exp_methods)

    n_matrices = len(find_common_matrices(best_dfs, exp_methods))
    print(f"Found {n_matrices} common matrices")

    counts = count_best_method(best_dfs, exp_methods, parameter = "VBR_nzblocks_count")
    print("BEST COUNT ON COMMON MATRICES: ", counts)
    make_plot(values_dict=counts, 
              title=f"Best method (lowest number of nonzero blocks) for on {n_matrices} matrices for {exp_methods}",
              ylabel="# of times X is the best method",
              save_path= f"{output_plot_dir}/blocks_best_count_plot_{exp_methods_string}")


    speedups = calculate_improvement(best_dfs, exp_methods[1:], parameter="VBR_nzblocks_count")
    print("SPEEDUPS GEOMETRIC MEAN ON COMMON MATRICES: ", speedups)
    make_improvements_barplot(best_dfs, 
                        exp_methods[1:],
                        parameter="VBR_nzblocks_count",
                        title=f"Median, 25th and 75th percentile of block number reduction on {n_matrices} matrices for {exp_methods}",
                        ylabel="Speed-up against non-reordered matrix",
                        save_path= f"{output_plot_dir}/blocks_speedup_median_{exp_methods_string}")

    make_improvements_violin(best_dfs, 
                        exp_methods[1:],
                        parameter="VBR_nzblocks_count",
                        title=f"Median, 25th and 75th percentile of block number reduction on {n_matrices} matrices for {exp_methods}",
                        ylabel="Block number reduction against non-reordered matrix",
                        save_path= f"{output_plot_dir}/block_speedup_violin_{exp_methods_string}")


    plot_improvements_distribution(best_dfs, 
                              exp_methods, 
                              parameter="VBR_nzblocks_count",
                              title=f"Cumulative distribution of block number reduction for on {n_matrices} matrices \n ({exp_methods})",

                              save_path= f"{output_plot_dir}/block_cumulative_hist_{exp_methods_string}.png"
                              )
    
    order_by= "clubs"
    plot_improvement_by_matrix( best_dfs, 
                            order_by,
                            exp_methods,
                            parameter="VBR_nzblocks_count", 
                            show="reduction_percent",
                            ylim=[0,400],
                            ylabel="Storage cost (compared to original)",
                            title=f"Storage cost by matrix for on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/block_matrix_order_by_{order_by}_{exp_methods_string}.png"
                              )
    
    order_by= "metis-edge-cut"
    plot_improvement_by_matrix( best_dfs, 
                            order_by,
                            exp_methods, 
                            parameter="VBR_nzblocks_count", 
                            show="reduction_percent",
                            ylim=[0,400],
                            ylabel="Storage cost (compared to original)",
                            title=f"Block number by matrix on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/block_matrix_order_by_{order_by}_{exp_methods_string}.png"
                              )
    
#ONLY CLUBS
#_____________________________________________________________________________

exp_methods =  ["original", "clubs"]
matrices = rectangular_matrices
order_by= "clubs"

plot_improvement_by_matrix( best_dfs, 
                        order_by,
                        exp_methods,
                        parameter="VBR_nzblocks_count", 
                        matrices=matrices,
                        title=f"Block reduction by matrix for block number reduction on {n_matrices} matrices \n ({exp_methods})",
                        save_path= f"{output_plot_dir}/block_matrices_rect_order_by_{order_by}_{exp_methods_string}.png"
                            )

matrices = square_matrices
plot_improvement_by_matrix( best_dfs, 
                        order_by,
                        exp_methods,
                        parameter="VBR_nzblocks_count", 
                        matrices=matrices,
                        title=f"Block reduction by matrix for block number reduction on {n_matrices} matrices \n ({exp_methods})",
                        save_path= f"{output_plot_dir}/block_matrices_square_order_by_{order_by}_{exp_methods_string}.png"
                            )

matrix_plot_dir = f"{output_plot_dir}/parameter_studies/"
os.makedirs(matrix_plot_dir, exist_ok=True)
method = "clubs"
for param in ["mask","centroid","tau"]:
    plot_club_performance_by_param(dfs_reordering, 
                                method=method,
                                param = param,
                                performance_param="VBR_nzblocks_count",
                                title=f'Performance variation under {param}', 
                                save_name=f"{matrix_plot_dir}/block_{method}_performance_variation_by_{param}")

#TRY: CALCULATE TOTAL INCLUDING RECTANGULAR FOR METIS










