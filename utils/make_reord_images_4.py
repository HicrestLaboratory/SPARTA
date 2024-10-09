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
parser.add_argument("--root_dir", nargs="?", type=str, default="results/results_02_10_2024/", help="the directory where the csv of the experiments are stored")
parser.add_argument("--bsize", nargs="?", type=str, default="64", help="size of the dense matrix")

args = parser.parse_args()

root_dir=args.root_dir
bsize=int(args.bsize)

if not os.path.isdir(root_dir):
    print(f"ERROR: Experiment directory {root_dir} does not exists.")
    exit(1)
#______________________________________________________


methods=["original", "clubs", "metis-edge-cut", "metis-volume","patoh"]
#methods=["original", "clubs", "metis-edge-cut","patoh"]

routines = ["spmmcsr", "spmmbsr","spmvcsr", "spmvbsr"]
#routines = ["spmmcsr", "spmmbsr"]


#subdirectories
csv_reordering_dir = root_dir + "reorder_csv"
csv_multiplication_dir = root_dir + "mult_csv"
output_plot_dir = root_dir + "reorder_plots"
os.makedirs(output_plot_dir, exist_ok=True)
for routine in routines:
    os.makedirs(f"{output_plot_dir}/{routine}", exist_ok=True)


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
    
    if "patoh" in method:
        dfs_reordering[method].rename(columns={'patoh_part': 'parts'}, inplace=True)

    dfs_reordering[method]['matrix'] = dfs_reordering[method]['matrix'].str.replace('_', '-', regex=False) #convention for matrix names
    
    #clean blocking data
    columns_to_drop = ["VBR_nzcount","VBR_average_height","VBR_longest_row"]
    dfs_reordering[method].drop(columns=columns_to_drop, inplace=True)
    dfs_reordering[method].rename(columns={"VBR_nzblocks_count" : "nnz_blocks"}, inplace=True)

#----------------------------------------------------------
#import reorder time data into dfs
#----------------------------------------------------------
for method in methods:
    reordering_file = f"{csv_reordering_dir}/clubs/reordering_time_results_clubs.csv"
    df_club_reordering_time = pd.read_csv(reordering_file, delim_whitespace=True, header=0)
    df_club_reordering_time['matrix'] = dfs_reordering[method]['matrix'].str.replace('_', '-', regex=False) #convention for matrix names
    df_club_reordering_time.rename(columns={"time" : "reordering_time"}, inplace=True)



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

    print(df_R.columns, dfs_reordering["original"].columns)

    #add original nz_blocks
    df_R_original = dfs_reordering["original"]
    df_R = df_R.merge(
                df_R_original[["matrix","nnz_blocks","rows","cols"]],
                on=["matrix","rows","cols"],
                how="left",
                suffixes=('', f'_original')
            )
    
    for routine in routines:

        #add original time
        df_M_original = dfs_mult[routine]["original"]
        df_R = df_R.merge(
                    df_M_original[["matrix","time"]],
                    on=["matrix"],
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






#______________________________________________________________________________________
# use patoh for failed CSR metis data
# Identify rows in df_metis where 'time_spmmcsr' is NaN
df_metis = dfs_reordering["metis-edge-cut"]
#df_patch = dfs_reordering["metis-volume"]
df_patch = dfs_reordering["patoh"]


# Merge df_metis with df_patoh on 'matrix' and 'parts'
df_merged = df_metis.merge(
    df_patch[['matrix', 'parts', 'time_spmmcsr']],
    on=['matrix', 'parts'],
    how='left',
    suffixes=('', '_patoh')
)

# Fill missing 'time_spmmcsr' in df_metis with values from df_patoh
df_merged['time_spmmcsr'] = df_merged['time_spmmcsr'].fillna(df_merged['time_spmmcsr_patoh'])

# Drop the extra 'time_spmmcsr_patoh' column
df_merged.drop(columns=['time_spmmcsr_patoh'], inplace=True)

# If needed, update df_metis with the merged data
dfs_reordering["metis-edge-cut"] = df_merged
#______________________________________________________________________________________




#calculate improvements and speedups
for method in methods:
    df_R = dfs_reordering[method]

    #calculate nz_blocks ratio
    df_R["blocks_ratio"] = df_R["nnz_blocks"]/df_R["nnz_blocks_original"]
    df_R["inverse_blocks_ratio"] = df_R["nnz_blocks_original"]/df_R["nnz_blocks"]
    df_R["density"] = df_R["nnz"]/(df_R["nnz_blocks"]*(64**2))

    for routine in routines:
        df_R[f"speedup_{routine}"] = df_R[f"time_{routine}_original"]/df_R[f"time_{routine}"]
        df_R[f"original_{routine}_failed"] = df_R[f"time_{routine}_original"].isna()
        df_R[f"method_{routine}_failed"] = (df_R[f"time_{routine}"].isna())
        df_R[f"gigaflops_{routine}"] = 10**(-8)*df_R["nnz"]*bsize/df_R[f"time_{routine}"]

        df_R.loc[df_R[f"original_{routine}_failed"] & ~df_R[f"method_{routine}_failed"], f"speedup_{routine}"] = float("inf")
        df_R.loc[~df_R[f"original_{routine}_failed"] & df_R[f"method_{routine}_failed"], f"speedup_{routine}"] = float("-inf")

        mult_successes = df_R[~df_R[f"time_{routine}"].isna()]["matrix"].unique()
        print(f"METHOD: {method}, ROUTINE: {routine}, SUCCESS MULT:{len(mult_successes)}")

print("ROWS IN PATOH:", method, dfs_reordering["patoh"].shape)


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

print("ROWS IN PATOH:", method, dfs_reordering["patoh"].shape)
print("ALL MATRICES: ", len(all_matrices_set))
print("SQUARE MATRICES: ", len(square_matrices_set))
print("RECTANGULAR MATRICES: ", len(rectangular_matrices_set))

print("RECTANGULAR MATRICES:", rectangular_matrices_set)
print("A rectangular matrix: ", rectangular_matrices_set.pop())
for routine in routines:
    print("COMMON MATRICES: ", routine, len(common_matrices_set[routine]))
    for method in methods:
        filtered_df = dfs_reordering[method][~dfs_reordering[method][f"method_{routine}_failed"]]
        unique_mats = filtered_df["matrix"].unique()
        unique_mats = dfs_reordering[method][~dfs_reordering[method][f"method_{routine}_failed"]]["matrix"].unique()
        unique_rect = len(set(unique_mats) & rectangular_matrices_set)
        unique_square = len(set(unique_mats) & square_matrices_set)
        print(f"***{method},{routine}: {len(unique_mats)} successes, of which RECT: {unique_rect} SQUARE: {unique_square}")
        #if "metis" in method:
            #print(f"MISSING METIS MATRICES FOR {routine}: {square_matrices_set - set(unique_mats)}")

#----------------------------------------------------------
# IMAGES
#----------------------------------------------------------


plt.rcParams.update({'font.size': 14})  # General font size
plt.xlabel('X Label', fontsize=16)     # X-axis label font size
plt.ylabel('Y Label', fontsize=16)     # Y-axis label font size
plt.title('Plot Title', fontsize=18, fontweight='bold')  # Title font size and bold
reordering_time_comparison(dfs_reordering["clubs"],df_club_reordering_time)

best_barplot(dfs_reordering, square_matrices_set, rectangular_matrices_set, methods, 
                 parameter = "nnz_blocks", 
                 ylabel = f"# of Matrices (Highest Block Density)",
                 fumbles=True,
                 fumbles_parameter = "inverse_blocks_ratio",
                 save_path=f"{output_plot_dir}/best_barplot_nnz_blocks_{bsize}")


for routine in routines:
    plot_parameter = "mask"
    method = "clubs"
    best_barplot_parameter(dfs_reordering,
                            square_matrices_set, 
                            rectangular_matrices_set, 
                            method = method, 
                            plot_parameter = plot_parameter, 
                            improvement_parameter = f"time_{routine}",
                            ratio_parameter=f"speedup_{routine}",
                            ylabel = "", 
                            xlabel= "Mask Size for CluB",
                            save_path=f"{output_plot_dir}/{routine}/{routine}_{plot_parameter}_time_best_plot_{method}")


plot_parameter = "mask"
method = "clubs"
#plot_params_values = [1,16,64]
plot_params_values = [1,16,64,256]
best_barplot_parameter(dfs_reordering,
                        square_matrices_set, 
                        rectangular_matrices_set, 
                        plot_params_values = plot_params_values,
                        method = method, 
                        plot_parameter = plot_parameter, 
                        improvement_parameter = f"nnz_blocks",
                        ratio_parameter="inverse_blocks_ratio",
                        ylabel = "", 
                        xlabel= "Mask Size for CluB",
                        save_path=f"{output_plot_dir}/{plot_parameter}_nnzb_best_plot_{method}_{bsize}")



for routine in routines:
    counts = count_best_method(dfs_reordering, square_matrices_set, parameter= f"time_{routine}" )
    print(f"BEST COUNT: {routine}", counts)
    best_barplot(dfs_reordering, square_matrices_set, rectangular_matrices_set, methods, 
                 parameter = f"time_{routine}", 
                 ylabel = "# of Matrices (Highest SpMM Speedup)",
                 fumbles = True,
                 fumbles_parameter = f"speedup_{routine}",
                 save_path=f"{output_plot_dir}/{routine}/{routine}_best_barplot_time")
    
for routine in routines:

    matrix_set = common_matrices_set[routine]
    make_improvements_barplot_and_distribution_2(dfs_reordering=dfs_reordering, 
                            methods=methods,
                            matrices=matrix_set,
                            ylabel=f"{routine} Speedup (Median)",
                            parameter=f"speedup_{routine}",
                            save_path=f"{output_plot_dir}/{routine}/{routine}_speedup_median_dist_time_common_matrices"
                            )
    

    matrix_set = square_matrices_set
    make_improvements_barplot_and_distribution_2(dfs_reordering=dfs_reordering, 
                            methods=methods,
                            matrices=matrix_set,
                            allow_missing=True,
                            ylabel=f"{routine_labels[routine]} Speedup (Median)",
                            parameter=f"speedup_{routine}",
                            save_path=f"{output_plot_dir}/{routine}/{routine}_speedup_median_dist_time_square_matrices"
                            )



for method in methods:
    compare_with="clubs"
    parameter = "blocks_ratio"
    ylabel="Relative BSR size"
    xlabel=f"Matrix ID (Sorted by {labels_dict[method]} Relative Size)"
    ylim=[0,4]
    yFormatter = percent_formatter
    plot_improvement_by_matrix(dfs_reordering,
                               methods= [compare_with,method],
                               order_by=method,
                               parameter=parameter,
                               ylabel=ylabel,
                               xlabel=xlabel,
                               ylim=ylim,
                               y_scale = "log",
                               original_line_y = 1,
                               yFormatter = yFormatter,
                               min_best=True,
                               matrices = square_matrices_set,
                               save_path=f"{output_plot_dir}/matrix_id_curve_{compare_with}by{method}_nnz_blocks_{bsize}")

    for routine in routines:
        compare_with="clubs"
        parameter = f"speedup_{routine}"
        ylabel="Speedup after reordering"
        xlabel=f"Matrix ID (sorted by {labels_dict[method]} speedup)"

        ylim=[0.75,2]
        if routine == "spmmbsr": ylim=[0.25,2]

        yFormatter = percent_improvement_formatter
        plot_improvement_by_matrix(dfs_reordering,
                                methods= [compare_with,method],
                                order_by=method,
                                parameter=parameter,
                                ylabel=ylabel,
                                xlabel=xlabel,
                                ylim=ylim,
                                yFormatter = yFormatter,
                                min_best=False,
                                matrices = square_matrices_set,
                                save_path=f"{output_plot_dir}/{routine}/{routine}_speedup_matrix_id_curve_{compare_with}by{method}")

        if False:
            speedup_vs_nnz_ratio(dfs_reordering,
                                method=method,  
                                x_parameter="nnz_blocks",
                                y_parameter=f"time_{routine}", 
                                matrices = square_matrices_set, 
                                xlim = [0,2], 
                                ylim=[0, 2], 
                                xlabel="Number of Nonzero Blocks (64 x 64)",
                                ylabel=f"{routine_labels[routine]} Time (ms)", 
                                save_path=f"{output_plot_dir}/{routine}/{routine}_time_vs_nnzblocks_{method}_{bsize}")

        speedup_vs_nnz_ratio(dfs_reordering,
                            method=method,  
                            x_parameter="density",
                            y_parameter=f"gigaflops_{routine}", 
                            matrices = square_matrices_set, 
                            xlim = [0,2], 
                            ylim=[0, 2],
                            xscale="linear",
                            yscale="log", 
                            xlabel="Density Within Nonzero Area",
                            ylabel=f"{routine_labels[routine]} GFLOPs", 
                            save_path=f"{output_plot_dir}/{routine}/{routine}_gflops_vs_density_{method}_{bsize}")

