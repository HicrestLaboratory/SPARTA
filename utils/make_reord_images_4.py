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


#----------------------------------------------------------
#import mult data into dfs
#----------------------------------------------------------
dfs = {}
for routine in routines:
    dfs[routine] = {method: pd.DataFrame() for method in methods}
    for method in methods:
        file_dir = f"{csv_multiplication_dir}/{routine}/{method}/"
        multiplication_file = file_dir + os.listdir(file_dir)[0]
        dfs[routine][method] = pd.read_csv(multiplication_file, delim_whitespace=True, header=0)
        df = dfs[routine][method]

        if routine == "spmmbsr":
            df.rename(columns={"nnz": "nnz_blocks"}, inplace = True)
            df.rename(columns={"rows": "rows_blocks"}, inplace = True)
            df.rename(columns={"cols": "cols_blocks"}, inplace = True)
            df["rows"] = df["rows_blocks"]
            df["cols"] = df["cols_blocks"]
            df["nnz"] = -1

        if "metis" in method:
            df = df[df['rows'] == df['cols']] 

        df['matrix'] = df['matrix'].str.replace('_', '-', regex=False) #convention for matrix names


        dfs[routine][method] = df

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for routine in routines:
    for method in methods:
        mats_mult = set(dfs[routine][method]["matrix"].unique())
        mats_reord = set(dfs_reordering[method]["matrix"].unique())
        mats_common = mats_mult & mats_reord
        only_mult = mats_mult - mats_common
        only_reord = mats_reord - mats_common

        #dfs[routine][method] = dfs[routine][method][~dfs[routine][method]['matrix'].isin(only_mult)]


        mats_mult = set(dfs[routine][method]["matrix"].unique())
        mats_reord = set(dfs_reordering[method]["matrix"].unique())
        mats_common = mats_mult & mats_reord
        only_mult = mats_mult - mats_common
        only_reord = mats_reord - mats_common

        print(f"{routine}, {method} : mult {len(mats_mult)}, reord {len(mats_reord)}, common {len(mats_common)}")
        print(sorted(list(only_mult))[:10])
        print(sorted(list(only_reord))[:10])
exit()


#Matrix df
all_matrices = pd.DataFrame()
for method in methods:
        df = dfs["spmmcsr"][method]
        all_matrices = pd.concat([all_matrices, df[['matrix', 'rows', 'cols','nnz']]])
        all_matrices = all_matrices.drop_duplicates()  # Ensure uniqueness

all_matrices["square"] = all_matrices["rows"] == all_matrices["cols"]
print(all_matrices)
exit()
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#add the speedups to the df
#----------------------------------------------------------
for routine in routines:
    for method in methods:
        df = dfs[routine][method]
        cols_to_merge = ["routine","matrix","time"]
        if routine == "spmmbsr": cols_to_merge.append("nnz_blocks")

        df = df.merge(
            dfs[routine]["original"][cols_to_merge],
            on=["routine","matrix"],
            how="outer",
            suffixes=('', f'_original')
        )

        df['speedup'] = df['time_original'] / df['time']
        df["original_failed"] = df["time_original"].isna()
        df["method_failed"] = (df["original_failed"] == False) & (df["speedup"].isna())

        if routine == "spmmbsr": df["block_improvement"] = df["nnz_blocks_original"] / df["nnz_blocks"]

        dfs[routine][method] = df



#add missing matrices
for routine in routines: 
    for method in methods:
        df = dfs[routine][method]
        df_csr = dfs["spmmcsr"]["original"]
        
        df = df.merge(df_csr[["matrix",]],
                            on = ["matrix"],
                            how = "outer",)
        df["original_failed"] = df["time_original"].isna()
        df["method_failed"] = df["time"].isna()        
        dfs[routine][method] = df



#add the row, cols, and nnz information to the bsr data
for method in methods:
    cols_to_merge = ["matrix","rows","cols","nnz"]
    df_bsr = dfs["spmmbsr"][method]

    df_csr = dfs["spmmcsr"]["original"]
    df_bsr = df_bsr.merge(df_csr[cols_to_merge],
                        on = ["matrix"],
                        how = "left",
                        suffixes=('', f'_csr'))
    df_bsr[["rows","cols","nnz"]] = df_bsr[["rows_csr","cols_csr","nnz_csr"]]
    df_bsr.drop(columns = ["rows_csr","cols_csr","nnz_csr"], inplace = True)
    dfs["spmmbsr"][method] = df_bsr

#identify square matrices
for routine in routines:
    for method in methods:
        df = dfs[routine][method]
        dfs[routine][method]["square"] = dfs[routine][method]["rows"] == dfs[routine][method]["cols"]




for routine in routines:
    for method in methods:
        df = dfs[routine][method]
        print(f"*****************************{routine}, {method}")
        square_mats = len(df[df["square"]]["matrix"].unique())
        rectangular_mats = len(df[df["square"] == False]["matrix"].unique())
        print("MATRICES: ", len(df["matrix"].unique()), " square: ", square_mats, " rect: ", rectangular_mats)
        print("METHOD FAILURES:", len(df[df["method_failed"]]["matrix"].unique()))
        print("ORIGINAL FAILURES:", len(df[df["original_failed"]]["matrix"].unique()))
        #print(dfs[routine][method])



#columns have been added to dfs[routine][method]: time_original and speedup
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
exit()




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










