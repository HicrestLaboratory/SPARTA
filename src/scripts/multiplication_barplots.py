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


bar_style = {}

bar_style["VBR-reord"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "orange"
}

bar_style["VBR-no-reord"] = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "blue"
}

bar_style["BELLPACK-no-reord"] = {
    "hatch" : "/",
    "edgecolor" : "black",
    "color" : "green"
}

bar_style["CSR"] = {
    "hatch" : "..",
    "edgecolor" : "black",
    "color" : "white"
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


def barplot(x_labels, x_ax_label, ys, y_labels, y_styles = {} , y_ax_label = "", savename = "test_barplot", title = ""):
    plt.figure()
    plt.xlabel(x_ax_label)
    plt.ylabel(y_ax_label)
    
    bars = len(ys)
    tot_width = 0.9
    barpos = -tot_width/2
    increment = tot_width/bars
    width = increment*0.95

    for y, label, style in zip(ys, y_labels, y_styles):
        print(len(y),label)
        x_pos = np.arange(barpos,len(x_labels) + barpos)
        plt.bar(x_pos,y,label=label, width = width, **style)
        barpos += increment

    plt.legend()
    plt.grid("both")
    plt.title(title)
    plt.xticks(range(len(x_labels)), x_labels, rotation=45)
    plt.savefig(savename + ".png",  bbox_inches='tight', dpi = 300)
    plt.close()    


exps = {}
exps["VBR-reord"] = (6,5)
exps["VBR-no-reord"] = (6,2)
exps["BELLPACK-no-reord"] = (3,2)
exps["CSR"] = (2,3)

def make_barplot(df, image_folder,B_cols, row_block_size, col_block_size):
    tmp_df = df.loc[(df["col_block_size"]==col_block_size) & (df["row_block_size"] == row_block_size) & (df["b_cols"] == B_cols)]
    matrices_names = [val.split("/")[-1].split(".")[0] for val in tmp_df["matrix"].unique()]

    data_lines = {exp_name : [] for exp_name in exps}
    if col_block_size != row_block_size: data_lines.pop("BELLPACK-no-reord")

    for matrix in tmp_df["matrix"].unique():
        for exp_name in data_lines.keys():
            M_algo, B_algo = exps[exp_name]
            values = tmp_df.loc[(tmp_df["matrix"] == matrix) & (tmp_df["multiplication_algo"] == M_algo) & (tmp_df["blocking_algo"] == B_algo)]["Speed-up against cuSparse"].values
            if len(values) == 0:
                data_lines[exp_name].append( 0 )
            else:
                data_lines[exp_name].append( values[0] )

    savename = f"{image_folder}/SpMM_time_barplot_{row_block_size}_{col_block_size}_{B_cols}"

    styles = [bar_style[exp] for exp in data_lines.keys()]

    barplot(
                x_labels = matrices_names,
                x_ax_label="graphs",
                ys=data_lines.values(),
                y_labels=data_lines.keys(),
                y_styles = styles,
                y_ax_label = "Speed-up against cuSparse", 
                savename = savename)

def make_barplot_best(df, image_folder,B_cols):
    tmp_df = df.loc[df["b_cols"] == B_cols]
    tmp_df = tmp_df.loc[tmp_df.groupby(["matrix","blocking_algo","multiplication_algo"])["avg_time_multiply"].idxmin()]
    matrices_names = [val.split("/")[-1].split(".")[0] for val in tmp_df["matrix"].unique()]

    data_lines = {exp_name : [] for exp_name in exps}

    for matrix in tmp_df["matrix"].unique():
        for exp_name in data_lines.keys():
            M_algo, B_algo = exps[exp_name]
            values = tmp_df.loc[(tmp_df["matrix"] == matrix) & (tmp_df["multiplication_algo"] == M_algo) & (tmp_df["blocking_algo"] == B_algo)]["Speed-up against cuSparse"].values
            if len(values) == 0:
                data_lines[exp_name].append( 0 )
            else:
                data_lines[exp_name].append( values[0] )

    savename = f"{image_folder}/SpMM_time_barplot_BEST_{B_cols}"

    styles = [bar_style[exp] for exp in data_lines.keys()]

    barplot(
                x_labels = matrices_names,
                x_ax_label="graphs",
                ys=data_lines.values(),
                y_labels=data_lines.keys(),
                y_styles = styles,
                y_ax_label = "Speed-up against cuSparse", 
                savename = savename)



def make_heatmap(df,image_folder,B_cols,exp_name):     
    M_algo, B_algo = exps[exp_name]
    colormap_variable = "Speed-up against cuSparse"
    heatmap_df = df[(df["multiplication_algo"]==M_algo) & (df["blocking_algo"]== B_algo) & (df["b_cols"] == B_cols)]
    table = heatmap_df.pivot_table(index="row_block_size", columns="col_block_size", values=colormap_variable, aggfunc='mean')
    table = table.sort_values(by=['row_block_size'], ascending=False)
    sns.heatmap(table,annot = True, cbar_kws={'label': exp_name})
    plt.ylabel("block height")
    plt.xlabel("block width")
    plt.savefig(f"{image_folder}/{exp_name}_heatmap_{colormap_variable}_{B_cols}.png",  bbox_inches='tight', dpi = 300)
    plt.close()

data_file = "test_suitsparse_3_multiplication.csv"
exp_name = "suitsparse_3_mult"
image_folder = "images/multiplication_images"
try: os.mkdir(image_folder) 
except: 1

df = pd.read_csv(data_file)
#df = get_dataframe_folder("results/suitsparse_collection_3")


df_CSR = df[df["multiplication_algo"] == 2][["matrix","b_cols","avg_time_multiply"]]
df = pd.merge(df,df_CSR, how = "left",on = ["matrix","b_cols"], suffixes=('','_CSR') )
df["Speed-up against cuSparse"] = df["avg_time_multiply_CSR"]/df["avg_time_multiply"]
#df_BELLPACK = df[df["multiplication_algo"] == 3][["matrix","b_cols","avg_time_multiply","row_block_size","col_block_size"]]
#df = pd.merge(df,df_BELLPACK, how = "left",on = ["matrix","b_cols","row_block_size","col_block_size"], suffixes=('','_BELLPACK') )
#df_VBR_no_reord = df.loc[(df["multiplication_algo"] == 6) & df["blocking_algo"] == 2][["matrix","b_cols","avg_time_multiply","row_block_size","col_block_size"]]
#df = pd.merge(df,df_BELLPACK, how = "left",on = ["matrix","b_cols","row_block_size","col_block_size"], suffixes=('','_VBR_no_reorder') )



df["density"] = df["nonzeros"].values/(df["rows"].values * df["cols"].values)
#df["block_density"] = df["nonzeros"].values/df["VBR_nzcount"].values
#df["dense_amp"] = df["block_density"].values/df["density"].values
df.sort_values(by=['density','matrix'], inplace=True)

#PREPARE HEATMAP FOR ALL BLOCK-SIZES


for B_cols in df["b_cols"].unique():
    make_barplot_best(df, image_folder,B_cols)
    for row_block_size in (32,64,128,512,1024):
        for col_block_size in (32,64,128,512,1024):
            make_barplot(df, image_folder,B_cols, row_block_size,col_block_size)

for B_cols in df["b_cols"].unique():
        for exp_name in exps:
            make_heatmap(df,image_folder,B_cols,exp_name)