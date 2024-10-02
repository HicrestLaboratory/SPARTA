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
parser.add_argument("--routine", nargs="?", type=str, default="spmmcsr", help="the multiplication routine to report on, e.g. spmmcsr")
parser.add_argument("--root_dir", nargs="?", type=str, default="results/results_02_10_2024/", help="the directory where the csv of the experiments are stored")

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


methods=["original", "clubs", "metis-edge-cut", "metis-volume", "patoh"]

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
degree_cutoff = 2
nnz_cutoff = 0
for method in methods: 
    #remove unorderable matrices
    dfs[method] = dfs[method].loc[dfs[method]["nnz"] >= degree_cutoff*dfs[method]["rows"]]
    dfs[method] = dfs[method].loc[dfs[method]["nnz"] >= nnz_cutoff]


all_matrices = set()
for method in methods:
    all_matrices = all_matrices.union(set(dfs[method]["matrix"].unique()))
print(f"TOTAL MATRIX COUNT: {len(all_matrices)}")


#find matrices for which original exist
matrices_with_original = dfs["original"]['matrix'].unique()
for method in methods:
    #dfs[method] = dfs[method][(dfs[method]["matrix"]).isin(matrices_with_original)]
    allowed_matrices = len(dfs[method]["matrix"].unique())
    print(f"Method: {method}; Found data for {allowed_matrices} matrices")

#find matrices common among ALL methods
common_matrices=set(matrices_with_original)
for method in methods:
    if method != "original":
        common_matrices &= set(dfs[method]['matrix'].unique())
print(f"Found {len(common_matrices)} common matrices")

#find rectangular and square matrices

rectangular_matrices = set()
square_matrices = set()
for method in methods: 
    rectangular_matrices |= set(dfs[method][dfs[method]['rows'] != dfs[method]['cols']]["matrix"].unique())
    square_matrices |= set(dfs[method][dfs[method]['rows'] == dfs[method]['cols']]["matrix"].unique())
print(f"RECTANGULAR: {len(rectangular_matrices)}, SQUARE: {len(square_matrices)}")


failed_matrices = {}
for method in methods:
    valid_matrices = set(all_matrices)
    working_matrices = set(dfs[method]["matrix"].unique()) 
    if "metis" in method:
        valid_matrices -= rectangular_matrices
        working_matrices -= rectangular_matrices
    failed_matrices[method] = valid_matrices - working_matrices
    print(f"METHOD: {method}, INVALID: {len(all_matrices - valid_matrices)}, FAILED: {len(failed_matrices[method])}, SUCCESS = {len(working_matrices)}")    
    

#find matrices that exist only for a method
method_matrices = {method : set(dfs[method]["matrix"].unique()) - common_matrices for method in methods}
best_dfs = {}
for method in methods:
    best_dfs[method] = find_best(dfs[method])

comparisons = [methods, ["original", "clubs", "metis-edge-cut", "patoh"], ["original", "clubs"], ["original", "metis-edge-cut", "metis-volume"]]


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

    counts = count_best_method(best_dfs, exp_methods)
    print("BEST COUNT ON COMMON MATRICES: ", counts)
    make_plot(values_dict=counts, 
              title=f"Best reordering method for {routine} on {n_matrices} matrices for {exp_methods}",
              ylabel="# of times X is the best method",
              save_path= f"{output_plot_dir}/{routine}_best_count_plot_{exp_methods_string}")


    speedups = calculate_improvement(best_dfs, exp_methods[1:])
    print("SPEEDUPS GEOMETRIC MEAN ON COMMON MATRICES: ", speedups)
    make_improvements_barplot(best_dfs, 
                        exp_methods[1:],
                        title=f"Median, 25th and 75th percentile of {routine} speedup on {n_matrices} matrices for {exp_methods}",
                        ylabel="Speed-up against non-reordered matrix",
                        save_path= f"{output_plot_dir}/{routine}_speedup_median_{exp_methods_string}")

    make_improvements_violin(best_dfs, 
                        exp_methods[1:],
                        title=f"Median, 25th and 75th percentile of {routine} speedup on {n_matrices} matrices for {exp_methods}",
                        ylabel="Speed-up against non-reordered matrix",
                        save_path= f"{output_plot_dir}/{routine}_speedup_violin_{exp_methods_string}")

    
    #make_plot(values_dict=speedups, 
    #          title=f"Geometric mean of {routine} speedup on {n_matrices} matrices for {exp_methods}",
    #          ylabel="Speed-up against non-reordered matrix",
    #          percent=True,
    #          save_path= f"{output_plot_dir}/{routine}_speedup_geomean_{exp_methods_string}")

    total_times = calculate_total_times(best_dfs, exp_methods)
    print("TOTAL TIMES ON COMMON MATRICES:", total_times)
    make_plot(values_dict=total_times, 
              title=f"Total time to run {routine} on {n_matrices} matrices for {exp_methods}",
              ylabel=f"Total {routine} time",
              save_path= f"{output_plot_dir}/{routine}_total_times_{exp_methods_string}")

    print("TOTAL TIME USING BEST: ", calculate_total_best_time(best_dfs, exp_methods))

    plot_improvements_distribution(best_dfs, 
                              exp_methods, 
                              title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",

                              save_path= f"{output_plot_dir}/{routine}_cumulative_hist_{exp_methods_string}.png"
                              )

    param = "rows"
    plot_speedup_vs_matrix_param(  best_dfs, 
                            exp_methods, 
                            param=param,
                            title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_speedup_vs_{param}_{exp_methods_string}.png"
                              )

    param = "density"
    plot_speedup_vs_matrix_param(  best_dfs, 
                            exp_methods, 
                            param=param,
                            title=f"Cumulative distribution of speedups for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_speedup_vs_{param}_{exp_methods_string}.png"
                              )
    
    order_by= "clubs"
    plot_improvement_by_matrix( best_dfs, 
                            order_by,
                            exp_methods, 
                            title=f"Speedup by matrix for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_matrix_order_by_{order_by}_{exp_methods_string}.png"
                              )
    
    order_by= "metis-edge-cut"
    plot_improvement_by_matrix( best_dfs, 
                            order_by,
                            exp_methods, 
                            title=f"Speedup by matrix for {routine} on {n_matrices} matrices \n ({exp_methods})",
                            save_path= f"{output_plot_dir}/{routine}_matrix_order_by_{order_by}_{exp_methods_string}.png"
                              )
    
#ONLY CLUBS
#_____________________________________________________________________________

exp_methods =  ["original", "clubs"]
matrices = rectangular_matrices
order_by= "clubs"

plot_improvement_by_matrix( best_dfs, 
                        order_by,
                        exp_methods,
                        matrices=matrices,
                        title=f"Speedup by matrix for {routine} on {n_matrices} matrices \n ({exp_methods})",
                        save_path= f"{output_plot_dir}/{routine}_matrices_rect_order_by_{order_by}_{exp_methods_string}.png"
                            )

matrices = square_matrices
plot_improvement_by_matrix( best_dfs, 
                        order_by,
                        exp_methods,
                        matrices=matrices,
                        title=f"Speedup by matrix for {routine} on {n_matrices} matrices \n ({exp_methods})",
                        save_path= f"{output_plot_dir}/{routine}_matrices_square_order_by_{order_by}_{exp_methods_string}.png"
                            )

matrix_plot_dir = f"{output_plot_dir}/parameter_studies/"
os.makedirs(matrix_plot_dir, exist_ok=True)
method = "clubs"
for param in ["mask","centroid","tau"]:
    plot_club_performance_by_param(dfs, 
                                method=method,
                                param = param, 
                                title=f'Performance variation under {param}', 
                                save_name=f"{matrix_plot_dir}/{routine}_{method}_performance_variation_by_{param}")

#TRY: CALCULATE TOTAL INCLUDING RECTANGULAR FOR METIS










