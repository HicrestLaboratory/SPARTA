import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
import argparse
from images_utils_2 import *

#----------------------ARGUMENTS
#______________________________________________________
parser = argparse.ArgumentParser(description="Analysis of multiplication times after reordering")
parser.add_argument("--root_dir", nargs="?", type=str, default="results/results_10-08-2024/", help="the directory where the csv of the experiments are stored")
parser.add_argument("--bsize", nargs="?", type=str, default="64", help="size of the dense matrix")

args = parser.parse_args()

root_dir=args.root_dir
bsize=int(args.bsize)

if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} does not exists.")
    exit(1)
#______________________________________________________


methods=["original", "clubs", "metis-edge-cut", "metis-volume"]
routines = ["spmmcsr", "spmmbsr"]

#subdirectories
csv_reordering_dir = root_dir + "reorder_csv"
csv_multiplication_dir = root_dir + "mult_csv"
output_plot_dir = root_dir + "reorder_plots/"
os.makedirs(output_plot_dir, exist_ok=True)


#----------------------------------------------------------
#import reorder data into dfs
#----------------------------------------------------------
dfs_reordering = {method: pd.DataFrame() for method in methods}
for method in methods:
    reordering_file = f"{csv_reordering_dir}/{method}_scramble0_bsize{bsize}.txt"
    dfs_reordering[method] = pd.read_csv(reordering_file, delim_whitespace=True, header=0)
    dfs_reordering[method].rename(columns={'matrix_name': 'matrix'}, inplace=True)
    if "metis" in method:
        dfs_reordering[method].rename(columns={'metis_part': 'parts'}, inplace=True)
        dfs_reordering[method].rename(columns={'metis_obj': 'objective'}, inplace=True)
        dfs_reordering[method] = dfs_reordering[method][dfs_reordering[method]['rows'] == dfs_reordering[method]['cols']] 
    dfs_reordering[method]['matrix'] = dfs_reordering[method]['matrix'].str.replace('_', '-', regex=False) #convention for matrix names
    dfs_reordering[method].rename(columns={'metis_part': 'parts'}, inplace=True)
    
    #clean blocking data
    columns_to_drop = ["VBR_nzcount","VBR_average_height","VBR_longest_row"]
    dfs_reordering[method].drop(columns=columns_to_drop, inplace=True)
    dfs_reordering[method].rename(columns={"VBR_nzblocks_count" : "nnz_blocks"}, inplace=True)


#----------------------------------------------------------
#import mult data into dfs
#----------------------------------------------------------
dfs_mult = {}
for routine in routines:
    dfs_mult[routine] = {method: pd.DataFrame() for method in methods}
    for method in methods:
        file_dir = f"{csv_multiplication_dir}/{routine}/{method}/"
        multiplication_file = file_dir + os.listdir(file_dir)[0]
        dfs_mult[routine][method] = pd.read_csv(multiplication_file, delim_whitespace=True, header=0)
        df = dfs_mult[routine][method]

        if "metis" in method:
            df = df.loc[df['rows'] == df['cols']].copy() 

        df.drop(columns = ["rows","cols","nnz"], inplace=True)

        df['matrix'] = df['matrix'].str.replace('_', '-', regex=False) #convention for matrix names

        dfs_mult[routine][method] = df
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#clean excess data (matrices not in suitesparse)
for routine in routines:
    for method in methods:
        mats_mult = set(dfs_mult[routine][method]["matrix"].unique())
        mats_reord = set(dfs_reordering[method]["matrix"].unique())
        mats_common = mats_mult & mats_reord
        only_mult = mats_mult - mats_common
        only_reord = mats_reord - mats_common

        dfs_mult[routine][method] = dfs_mult[routine][method][~dfs_mult[routine][method]['matrix'].isin(only_mult)]

        mats_mult = set(dfs_mult[routine][method]["matrix"].unique())
        mats_reord = set(dfs_reordering[method]["matrix"].unique())
        mats_common = mats_mult & mats_reord
        only_mult = mats_mult - mats_common
        only_reord = mats_reord - mats_common

        print(f"{routine}, {method} : mult {len(mats_mult)}, reord {len(mats_reord)}, common {len(mats_common)}")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#add mult times to reordering
for method in methods:
    df_R = dfs_reordering[method]

    #add original nz_blocks
    df_R_original = dfs_reordering["original"]
    df_R = df_R.merge(
                df_R_original[["matrix","nnz_blocks"]],
                on=["matrix",],
                how="left",
                suffixes=('', f'_original')
            )

    for routine in routines:

        #add original time
        df_M_original = dfs_mult[routine]["original"]
        df_R = df_R.merge(
                    df_M_original[["matrix","time"]],
                    on=["matrix",],
                    how="left",
                )
        df_R.rename(columns={f'time': f'time_{routine}_original'}, inplace=True)

        #add routine results
        df_M = dfs_mult[routine][method]
        common_columns = list(set(df_M.columns) & set(df_R.columns)) 
        df_R = df_R.merge(
                    df_M[common_columns + ["time",]],
                    on=common_columns,
                    how="left",
                )
        df_R.rename(columns={f'time': f'time_{routine}'}, inplace=True)
        
    dfs_reordering[method] = df_R


#calculate improvements and speedups
for method in methods:
    df_R = dfs_reordering[method]

    #calculate nz_blocks ratio
    df_R["blocks_ratio"] = df_R["nnz_blocks"]/df_R["nnz_blocks_original"]
    for routine in routines:
        df_R[f"speedup_{routine}"] = df_R[f"time_{routine}_original"]/df_R[f"time_{routine}"]
        df_R[f"original_{routine}_failed"] = df_R[f"time_{routine}_original"].isna()
        df_R[f"method_{routine}_failed"] = (df_R[f"time_{routine}"].isna())

        mult_successes = df_R[~df_R[f"time_{routine}"].isna()]["matrix"].unique()
        print(f"METHOD: {method}, ROUTINE: {routine}, SUCCESS MULT:{len(mult_successes)}")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Remove matrices that fail METIS reordering.
df_metis = dfs_reordering["metis-edge-cut"]
successes_in_metis_reordering = set(df_metis["matrix"].unique())
failed_matrices = set(dfs_reordering["original"][dfs_reordering["original"]["rows"] ==  dfs_reordering["original"]["cols"]]["matrix"]) - successes_in_metis_reordering
print("FAILURES IN METIS REORDERING:", len(failed_matrices))

for method in methods:
    dfs_reordering[method] = dfs_reordering[method].loc[~dfs_reordering[method]['matrix'].isin(failed_matrices)]
    print(f"removed METIS-FAILED matrices from METHOD {method}, now at {len(dfs_reordering[method]['matrix'].unique())}")


#----------------------------------------------------------#
# MATRIX SETS
#----------------------------------------------------------
all_matrices = dfs_reordering["original"].loc[:,["matrix","rows","cols","nnz"]]
all_matrices["square"] = all_matrices["rows"] == all_matrices["cols"]
all_matrices_set = set(all_matrices["matrix"].unique())
square_matrices_set = set(all_matrices[all_matrices["square"]]["matrix"].unique())
rectangular_matrices_set = set(all_matrices[~all_matrices["square"]]["matrix"].unique())
common_matrices_set = {routine: set(dfs_reordering["original"]["matrix"]) for routine in routines}
for routine in routines:
    for method in methods:
        df = dfs_reordering[method][~dfs_reordering[method][f"time_{routine}"].isna()]
        matrices_in_df = set(df["matrix"])
        common_matrices_set[routine] &= matrices_in_df

print("ALL MATRICES: ", len(all_matrices_set))
print("SQUARE MATRICES: ", len(square_matrices_set))
print("RECTANGULAR MATRICES: ", len(rectangular_matrices_set))
print("A rectangular matrix: ", rectangular_matrices_set.pop())
for routine in routines:
    print("COMMON MATRICES: ", routine, len(common_matrices_set[routine]))

#----------------------------------------------------------
# IMAGES
#----------------------------------------------------------

best_barplot(dfs_reordering, square_matrices_set, rectangular_matrices_set, methods, 
                 parameter = "nnz_blocks", 
                 ylabel = "# of times reordering results in fewest nonzero blocks",
                 save_path=f"{output_plot_dir}/best_barplot_{routine}_nnz_blocks_{methods}")

for routine in routines:
    counts = count_best_method(dfs_reordering, square_matrices_set, methods, parameter= f"time_{routine}" )
    print(f"BEST COUNT: {routine}", counts)
    best_barplot(dfs_reordering, square_matrices_set, rectangular_matrices_set, methods, 
                 parameter = f"time_{routine}", 
                 ylabel = "# of times reordering results in fastest multiplication",
                 save_path=f"{output_plot_dir}/best_barplot_{routine}_time_{methods}")
    

exit()


for exp_methods in comparisons:
    exp_methods_string = "_".join(exp_methods)
    print("*****************")
    print("Comparing: ", exp_methods)

    n_matrices = len(find_common_matrices(best_dfs, exp_methods))
    print(f"Found {n_matrices} common matrices")

    counts = count_best_method(best_dfs, exp_methods, parameter = "VBR_nzblocks_count")
    print("BEST COUNT ON COMMON MATRICES: ", counts)
    make_barplot(values_dict=counts, 
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










