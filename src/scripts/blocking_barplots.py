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
from matplotlib import cm, colors
import argparse



global_label_dict = {
    "block_density" : "Blocked density (IC)",
    "density" : "Original density",
    "nonzeros" : "# of nonzeros in the original matrix",
    "dense-amp" : "Density amplification (against unblocked)",
    "relative-dense-amp" : "Density amplification (against natural blocking)",
    "block ratio": "Shape (height / width)",
    "block_density_no_reord" : "Blocked density (natural)"
}

global_exp_dict = {}

global_exp_dict["reordered"] = {
    "style" : {
        "hatch" : "//",
        "edgecolor" : "black",
        "color" : "orange",
        "label" : "IC blocking"
    }
}

global_exp_dict["no-reordered"] = {
    "style": {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "blue",
    "label" : "natural blocking"
    }
}

global_exp_dict["original"] = {
    "style" : {
    "hatch" : "..",
    "edgecolor" : "black",
    "color" : "white",
    "label" : "orginal matrix"
    }
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

def barplot(x_labels, x_ax_label, ys, y_styles, y_ax_label, yscale = "log", savename = "test_barplot"):
    plt.figure()
    plt.xlabel(x_ax_label)
    plt.ylabel(y_ax_label)
    
    bars = len(ys)
    tot_width = 0.9
    barpos = -tot_width/2
    increment = tot_width/bars
    width = increment*0.95

    for y, style in zip(ys, y_styles):
        x_pos = np.arange(barpos,len(x_labels) + barpos)
        plt.bar(x_pos,y, width = width, **style)
        barpos += increment

    plt.legend()
    plt.yscale(yscale)
    plt.grid("both")
    plt.title(f"block height = {row_block_size}, block width = {col_block_size}")
    plt.xticks(range(len(x_labels)), x_labels, rotation=90)
    plt.savefig(savename + ".png",  bbox_inches='tight', dpi = 300)
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
                    default='images/blocking_images/test_images',
                    help='The directory where to save the images')
args = vars(parser.parse_args())


data_file = args["file"]
image_folder = args["imagedir"]
try: os.mkdir(image_folder) 
except: 1

df = pd.read_csv(data_file)
#df = get_dataframe_folder("results/suitsparse_collection_3")

d = (0.00005, 0.1)
n = (20000,100000)

df = df.loc[df.groupby(["matrix","blocking_algo","col_block_size","row_block_size"])["VBR_nzblocks_count"].idxmin()]

df["density"] = df["nonzeros"].values/(df["rows"].values * df["cols"].values)
df = df[(df["density"] >= d[0]) & (df["density"] <= d[1])]
df = df[(df["rows"] >= n[0]) & (df["rows"] <= n[1])]
df = df[(df["cols"] >= n[0]) & (df["cols"] <= n[1])]

df["block_density"] = df["nonzeros"].values/df["VBR_nzcount"].values
df["dense-amp"] = df["block_density"].values/df["density"].values
df["block_area"] = df["row_block_size"].values*df["col_block_size"].values
df["block ratio"] = df["row_block_size"].values/df["col_block_size"].values


df_VBR_no_reord = df[df["blocking_algo"] == 2][["matrix","col_block_size","row_block_size","block_density"]]
df = pd.merge(df,df_VBR_no_reord, how = "left",on = ["matrix","row_block_size","col_block_size"], suffixes=('','_no_reord') )
df["relative-dense-amp"] = df["block_density"]/df["block_density_no_reord"]
df.loc[df["relative-dense-amp"] < 1, "relative-dense-amp"] = 1

for block_size in (64,128,256,512,1024):
    print(f"BLOCK SIZE: {block_size}")
    print("BAD", len(df.loc[(df["relative-dense-amp"] == 1) & (df["row_block_size"] == block_size)]))
    print("GOOD", len(df.loc[(df["relative-dense-amp"] != 1) & (df["row_block_size"] == block_size)]))


df.sort_values(by=['density','matrix'], inplace=True)

np.set_printoptions(suppress = True)

print(df["matrix"].unique())
for row_block_size in (32,64,128,256,512,1024):
    for col_block_size in (64,128,256,512,1024):
        print(f"row_block_size, col_block_size = {row_block_size},{col_block_size}")
        tmp_df = df.loc[(df["blocking_algo"] == 5) & (df["col_block_size"] == col_block_size) & (df["row_block_size"] == row_block_size)]
        print(tmp_df["tau"].values)


def make_scatter(x_var = "density", y_var = "relative-dense-amp", row_block_size = 128, col_block_size = 128, xscale = "log", yscale = "linear", drawline = 0):
    blocking_algo = 5
    fig,ax = plt.subplots()
    tmp_df = df.loc[(df["blocking_algo"] == blocking_algo) & (df["row_block_size"] == row_block_size) & (df["col_block_size"] == col_block_size)]
    ax = tmp_df.plot.scatter(
                    x= x_var,
                    y= y_var,
                    colormap='viridis',alpha=0.5, edgecolor = "black")
    plt.ylabel(global_label_dict[y_var])
    if drawline:
        plt.axhline(drawline,color = "red", linestyle= "--")
    #plt.xlim(0.00001,1)
    #plt.ylim(0.00001,1)
    plt.xscale(xscale)
    plt.yscale(yscale)
    #plt.yscale("log")
    plt.xlabel(global_label_dict[x_var])
    plt.savefig(f"{image_folder}/scatterplot_{x_var}_vs_{y_var}_b{row_block_size}_B{col_block_size}_algo_{blocking_algo}.png", bbox_inches='tight', dpi = 300)
    plt.close()

for block_size in (64,128,256,512,1024):
        #try:
        make_scatter(x_var ="block_density_no_reord", y_var = "relative-dense-amp", row_block_size = block_size, col_block_size = block_size, drawline = 1)
        #make_scatter(x_var ="block_density_no_reord", y_var = "block_density", row_block_size = block_size, col_block_size = block_size)

        #except:
        #    print(f"COULD NOT MAKE SCATTER FOR {block_size} x {block_size}")


"""
x_var = "density"
y_var = "relative-dense-amp"
for row_block_size in (32,64,128,256,512,1024):
    for col_block_size in (32,64,128,256,512,1024):
        try:
            fig,ax = plt.subplots()
            tmp_df = df.loc[(df["blocking_algo"] == 5) & (df["row_block_size"] == row_block_size) & (df["col_block_size"] == row_block_size)]
            tmp_df.plot.scatter(
                            x= x_var,
                            y= y_var,
                            colormap='viridis',alpha=0.5, edgecolor = "black", marker = "s", ax = ax, label = "Blocking after reordering")
            plt.ylabel(global_label_dict[y_var])
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xscale("log")
            plt.savefig(f"{image_folder}/scatterplot_{x_var}_vs_{y_var}_{row_block_size}x{col_block_size}.png", bbox_inches='tight', dpi = 300)
        except:
            print(f"FAILED TO PRODUCE SCATTER PLOT FOR {row_block_size} x {col_block_size}")
        plt.close()
"""


color_var = "relative-dense-amp"
color_label = "Density Amplification (against natural blocking)"

plt.figure()
vmin = 1
vmax = 2
table = df.loc[(df["blocking_algo"] == 5) & (df["relative-dense-amp"] > 1)]
table = table.pivot_table(index="row_block_size", columns="col_block_size", values=color_var, aggfunc='mean')
table = table.sort_values(by = ["row_block_size",], ascending = False)
sns.heatmap(table,annot = True, cbar_kws={'label': color_label},vmin = vmin, vmax = vmax)
plt.ylabel("block height")
plt.xlabel("block width")
plt.savefig(f"{image_folder}/heatmap_{color_var}.png",  bbox_inches='tight', dpi = 300)
plt.close()
    

"""
for algo in (2,5):
    heatmap_df = df[df["blocking_algo"]==algo].pivot_table(index="row_block_size", columns="col_block_size", values=colormap_variable, aggfunc='mean')
    heatmap_df = heatmap_df.sort_values(by=['row_block_size'], ascending=False)
    sns.heatmap(heatmap_df,annot = True, cbar_kws={'label': color_label},vmin = 1, vmax = max_var_value)
    plt.ylabel("block height")
    plt.xlabel("block width")
    plt.savefig(f"{image_folder}/heatmap_{colormap_variable}_algo_{algo}.png",  bbox_inches='tight', dpi = 300)
    plt.close()
"""




for row_block_size in sorted(df["row_block_size"].unique()):
    for col_block_size in sorted(df["col_block_size"].unique()):
        try:
            tmp_df = df.loc[(df["col_block_size"]==col_block_size) & (df["row_block_size"] == row_block_size)]
            matrices_names = [val.split("/")[-1].split(".")[0] for val in tmp_df["matrix"].unique()]

            reorder_df = tmp_df[tmp_df["blocking_algo"] == 5]
            no_reorder_df = tmp_df[tmp_df["blocking_algo"] == 2]
            try:
                assert np.equal(reorder_df["matrix"].values, no_reorder_df["matrix"].values).all() #check that there are the same n of values for all matrices
                assert len(reorder_df) > 0
            except:
                print(f"MISSING DATA for blocks {row_block_size} x {col_block_size}")
                continue

            print(f"FOUND DATA for {row_block_size} x {col_block_size}")
            data_lines = {}

            data_lines["blocking, no reordering"] = no_reorder_df["block_density"].values
            data_lines["blocking, reordering"] = reorder_df["block_density"].values
            data_lines["no blocking (original)"] = tmp_df["density"].unique()   

            savename = f"{image_folder}/Density_barplot_{row_block_size}_{col_block_size}"

            exps = ""
            styles = [exp["style"] for exp in global_exp_dict]

            barplot(
                    x_labels = matrices_names,
                    x_ax_label="graphs",
                    ys=data_lines.values(),
                    y_labels=data_lines.keys(),
                    y_ax_label = "Density", 
                    savename = savename)
        except:
            print(f"NO BARPLOT FOR {row_block_size} x {col_block_size}")
    plt.close()
