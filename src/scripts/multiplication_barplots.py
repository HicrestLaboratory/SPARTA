# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:35:55 2022

@author: Paolo
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import operator
import pandas as pd
import os as os
import seaborn as sns
import argparse



global_label_dict = {
    "block_density" : "Blocked density",
    "density" : "Original density",
    "dense-amp" : "Density amplification vs natural blocking",
    "relative_dense_amp" : "Density amplification (against natural blocking)",
    "block ratio": "Shape (height / width)",
    "speed-vs-cusparse": "Speed-up against cuSparse-CSR",
    "speed-vs-no-reord": "Speed-up against cuBLAS-BCSR with natural blocking"
}

global_exp_dict = {}

global_exp_dict["VBR-reord"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "orange",
    "label" : "cuBLAS-BCSR (IC blocking)"

}

global_exp_dict["VBR-no-reord"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "blue",
    "label" : "cuBLAS-BCSR (natural blocking)"
}

global_exp_dict["BELLPACK-no-reord"] = {
    "hatch" : "/",
    "edgecolor" : "black",
    "color" : "green",
    "label" : "cuSparse-BELLPACK (natural blocking)"

}

global_exp_dict["BELLPACK-reord"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "green",
    "label" : "cuSparse-BELLPACK (IC blocking)"

}


global_exp_dict["GEMM"] = {
    "hatch" : "--",
    "edgecolor" : "black",
    "color" : "red",
    "label" : "cuBLAS-gemm"

}

global_exp_dict["CSR"] = {
    "hatch" : "..",
    "edgecolor" : "black",
    "color" : "white",
    "label" : "cuSparse-CSR"

}

global_exp_dict["CUTLASS_BELLPACK"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "gray",
    "label" : "CUTLASS-BELLPACK"

}


def get_dataframe_folder(folder):
    df = pd.DataFrame()
    for matrix_folder in glob.glob(f"{folder}/*"):
        for experiment in glob.glob(f"{matrix_folder}/*"):
            try:
                df1 = pd.read_csv(experiment, skiprows = 0)
                df = pd.concat([df, df1], ignore_index=True, sort=False)
            except:
                print(f"CANNOT LOAD {experiment}")
    if df.empty:
        print("ALTERNATIVE LOAD", folder)
        for experiment in glob.glob(f"{folder}/*"):
            print(experiment)
            try:
                df1 = pd.read_csv(experiment, skiprows = 0)
                df = pd.concat([df, df1], ignore_index=True, sort=False)
            except:
                print(f"CANNOT LOAD {experiment}")
    return df        

def get_line(df, constraints):
    query = ""
    for key in constraints:
        query += f"{key} == {constraints[key]} &"
    query = query[:-2]
    return df.query(query)


def barplot(x_labels, x_ax_label, ys, y_styles = {} , y_ax_label = "", title = ""):
    plt.xlabel(x_ax_label)
    plt.ylabel(y_ax_label)
    
    bars = len(ys)
    tot_width = 0.9
    increment = tot_width/bars
    barpos = -tot_width/2
    width = increment*0.95

    for y, style in zip(ys, y_styles):
        x_pos = np.arange(barpos,len(x_labels) + barpos)
        plt.bar(x_pos,y, width = width, **style)
        barpos += increment

    plt.grid("both")
    plt.title(title)
    plt.xticks(range(len(x_labels)), x_labels, rotation=90)
 


exps = {}
exps["VBR-no-reord"] = (6,2)
exps["VBR-reord"] = (6,5)
exps["BELLPACK-no-reord"] = (3,2)
#exps["BELLPACK-reord"] = (3,5)
exps["CSR"] = (2,3)
#exps["GEMM"] = (1,3)
#exps["CUTLASS_BELLPACK"] = (8,5)



def make_barplot(df, image_folder,B_cols, row_block_size, col_block_size, variable = "speed-vs-cusparse"):
    tmp_df = df.loc[(df["col_block_size"]==col_block_size) & (df["row_block_size"] == row_block_size) & (df["b_cols"] == B_cols)]
    matrices_names = [val.split("/")[-1].split(".")[0] for val in tmp_df["matrix"].unique()]

    data_lines = {exp_name : [] for exp_name in exps}
    data_lines.pop("CSR")
    if col_block_size != row_block_size: data_lines.pop("BELLPACK-no-reord")

    for matrix in tmp_df["matrix"].unique():
        for exp_name in data_lines.keys():
            M_algo, B_algo = exps[exp_name]
            values = tmp_df.loc[(tmp_df["matrix"] == matrix) & (tmp_df["multiplication_algo"] == M_algo) & (tmp_df["blocking_algo"] == B_algo)][variable].values
            if len(values) == 0:
                data_lines[exp_name].append( 0 )
            else:
                data_lines[exp_name].append( values[0] )

    savename = f"{image_folder}/SpMM_barplot_{variable}_{row_block_size}_{col_block_size}_{B_cols}"

    styles = [global_exp_dict[exp] for exp in data_lines.keys()]


    plt.figure()
    barplot(
                x_labels = matrices_names,
                x_ax_label="graphs",
                ys=data_lines.values(),
                y_styles = styles,
                y_ax_label = global_label_dict[variable]
                )
    plt.legend()
    plt.savefig(savename + ".png",  bbox_inches='tight', dpi = 300)
    plt.close()   
    

def make_barplot_best(df, image_folder,B_cols, variable = "speed-vs-cusparse"):
    tmp_df = df.loc[df["b_cols"] == B_cols]
    tmp_df = tmp_df.loc[df["row_block_size"] == df["col_block_size"]]
    tmp_df = tmp_df.loc[tmp_df.groupby(["matrix","blocking_algo","multiplication_algo"])["avg_time_multiply"].idxmin()]
    tmp_df.sort_values(by=['density','matrix'], inplace=True)
    matrices_names = [val.split("/")[-1].split(".")[0] for val in tmp_df["matrix"].unique()]

    data_lines = {exp_name : [] for exp_name in exps}
    data_lines.pop("CSR")

    for matrix in tmp_df["matrix"].unique():
        for exp_name in data_lines.keys():
            M_algo, B_algo = exps[exp_name]
            values = tmp_df.loc[(tmp_df["matrix"] == matrix) & (tmp_df["multiplication_algo"] == M_algo) & (tmp_df["blocking_algo"] == B_algo)][variable].values
            if len(values) == 0:
                data_lines[exp_name].append( 0 )
            else:
                data_lines[exp_name].append( values[0] )

    savename = f"{image_folder}/SpMM_time_barplot_BEST_{B_cols}"

    styles = [global_exp_dict[exp] for exp in data_lines.keys()]


    plt.figure()
    barplot(
                x_labels = matrices_names,
                x_ax_label="Sparse matrices",
                ys=data_lines.values(),
                y_styles = styles,
                y_ax_label = global_label_dict[variable])
    plt.axhline(linewidth=2, y = 1, color = "red", label = global_exp_dict["CSR"]["label"])
    plt.legend()
    plt.savefig(savename + ".png",  bbox_inches='tight', dpi = 300)
    plt.close()   





def make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable = "speed-vs-cusparse"):     
    M_algo, B_algo = exps[exp_name]
    heatmap_df = df[(df["multiplication_algo"]==M_algo) & (df["blocking_algo"]== B_algo) & (df["b_cols"] == B_cols)]
    table = heatmap_df.pivot_table(index="row_block_size", columns="col_block_size", values=colormap_variable, aggfunc='mean')
    table = table.sort_values(by=['row_block_size'], ascending=False)
    if colormap_variable in global_label_dict:
        c_label = global_label_dict[colormap_variable]
    else:
        c_label = colormap_variable
    sns.heatmap(table,annot = True, cbar_kws={'label': c_label})
    plt.ylabel("block height")
    plt.xlabel("block width")
    plt.savefig(f"{image_folder}/{exp_name}_heatmap_{colormap_variable}_{B_cols}.png",  bbox_inches='tight', dpi = 300)
    plt.close()


def make_scatter(df,B_cols,var_x = "dense-amp", var_y = "speed-vs-no-reord"):
    plt.figure()
    tmp_df = df.loc[(df["b_cols"] == B_cols) & (df["multiplication_algo"] == 6) & (df["blocking_algo"] == 5)]
    sns.scatterplot(data=tmp_df, x=var_x, y=var_y)
    plt.xlabel(global_label_dict[var_x])
    plt.ylabel(global_label_dict[var_y])
    plt.savefig(f"{image_folder}/scatter_plot_{var_x}_vs_{var_y}",  bbox_inches='tight', dpi = 300)
    plt.close()


parser = argparse.ArgumentParser(description='Plots for multiplication experiments')
parser.add_argument('-f',
                    '--file',
                    type=str,
                    default='results/test_suitsparse_3_multiplication.csv',
                    help='The csv with the results to be plotted')
parser.add_argument('-i',
                    '--imagedir',
                    type = str,
                    default='images/multiplication_images/',
                    help='The directory where to save the images')
args = vars(parser.parse_args())


data_file = args["file"]
image_folder = args["imagedir"]
try: os.mkdir(image_folder) 
except: 1

df = pd.read_csv(data_file)
#df = get_dataframe_folder("results/suitsparse_collection_3")

df_CSR = df[df["multiplication_algo"] == 2][["matrix","b_cols","avg_time_multiply"]]
df = pd.merge(df,df_CSR, how = "left",on = ["matrix","b_cols"], suffixes=('','_CSR') )
print(list(df["nonzeros"].unique()))
df["block_density"] = np.where(df["nonzeros"].str.isnumeric(), df["nonzeros"]/df["VBR_nzcount"],np.nan)

df_VBR_no_reord = df[(df["multiplication_algo"] == 6) & (df["blocking_algo"] == 2)][["matrix","b_cols","avg_time_multiply","col_block_size","row_block_size","block_density"]]

df = pd.merge(df,df_VBR_no_reord, how = "left",on = ["matrix","b_cols","row_block_size","col_block_size"], suffixes=('','_VBR_no_reord') )

df["dense-amp"] = df["block_density"]/df["block_density_VBR_no_reord"]
df["speed-vs-no-reord"] = df["avg_time_multiply_VBR_no_reord"]/df["avg_time_multiply"]
df["speed-vs-cusparse"] = df["avg_time_multiply_CSR"]/df["avg_time_multiply"]
df["time-per-block"] = df["avg_time_multiply"]/df["VBR_nzblocks_count"]
df["time-per-area"] = df["avg_time_multiply"]/df["VBR_nzcount"]
df["time-per-true-nonzero"] = df["avg_time_multiply"]/df["nonzeros"]

#df_BELLPACK = df[df["multiplication_algo"] == 3][["matrix","b_cols","avg_time_multiply","row_block_size","col_block_size"]]
#df = pd.merge(df,df_BELLPACK, how = "left",on = ["matrix","b_cols","row_block_size","col_block_size"], suffixes=('','_BELLPACK') )
#df_VBR_no_reord = df.loc[(df["multiplication_algo"] == 6) & df["blocking_algo"] == 2][["matrix","b_cols","avg_time_multiply","row_block_size","col_block_size"]]
#df = pd.merge(df,df_BELLPACK, how = "left",on = ["matrix","b_cols","row_block_size","col_block_size"], suffixes=('','_VBR_no_reorder') )


df["density"] = df["nonzeros"].values/(df["rows"].values * df["cols"].values)
#df["block_density"] = df["nonzeros"].values/df["VBR_nzcount"].values
#df["dense_amp"] = df["block_density"].values/df["density"].values
df.sort_values(by=['density','matrix'], inplace=True)

#PREPARE HEATMAP FOR ALL BLOCK-SIZES


print("Making barplots")
for B_cols in df["b_cols"].unique():
    make_scatter(df,B_cols)
    print(f"***** for B_cols = {B_cols}")
    make_barplot_best(df, image_folder,B_cols)
    for row_block_size in (32,64,128,256,512,1024):
        for col_block_size in (32,64,128,256,512,1024):
            make_barplot(df, image_folder,B_cols, row_block_size,col_block_size)
            1

print("Making heatmaps")
for B_cols in df["b_cols"].unique():
        for exp_name in exps:
            try:
                make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable= "time-per-true-nonzero")
                make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable= "speed-vs-cusparse")
                make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable= "speed-vs-no-reord")
                make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable= "time-per-block")
                make_heatmap(df,image_folder,B_cols,exp_name, colormap_variable= "time-per-area")
            except: 
                print(f"no heatmap for {exp_name}")