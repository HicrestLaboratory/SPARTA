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



label_dict = {
    "block_density" : "Blocked density",
    "density" : "Original density",
    "dense_amp" : "Density amplification (against unblocked)",
    "block ratio": "Shape (height / width)"
}


bar_style_reordering = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "blue"
}

bar_style_blocking = {
    "hatch" : "//",
    "edgecolor" : "black",
    "color" : "orange"
}

bar_style_natural = {
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





"""
folder_name = "suitsparse_collection_2"
folder = f"results/{folder_name}"


variable = "block_density"
variable_2 = "tau"
ylabel = "Density"
savename = f"results/images/suitsparse_collection/suitsparse_collection_{variable}"
df = pd.DataFrame()
for fold_num in (3,4):
    folder_name = "suitsparse_collection_" + str(fold_num)
    folder = f"results/{folder_name}"
    df = pd.concat([df,get_dataframe(folder)],ignore_index=True, sort=False)

constraints = {}
for row_block_size in sorted(df["row_block_size"].unique()):
    for col_block_size in sorted(df["col_block_size"].unique()):
        plt.figure()
        plt.xlabel("graphs")
        plt.ylabel(ylabel)
        
        var_values = {};
        bars = 3
        barpos = -0.45
        increment = 0.9/bars
        width = increment*0.95

        for algo,algoname in zip((2,5),("blocking, no-reordering","blocking, reordering")):
            constraints["col_block_size"] = col_block_size
            constraints["row_block_size"] = row_block_size
            constraints["blocking_algo"] = algo
            res_df = get_best_blockings(df, variable = "VBR_nzblocks_count", variable_2 = "tau", constraints = constraints)
            matrices = [val.split("/")[-1].split(".")[0] for val in res_df["matrix"].values]
            taus = res_df["tau"].values
            res_df["density"] = res_df["nonzeros"].values/(res_df["rows"].values * res_df["cols"].values)
            res_df["block_density"] = res_df["nonzeros"].values/res_df["VBR_nzcount"].values
            var_values[algo] = res_df[variable].values
            print(algo, row_block_size, col_block_size, matrices, taus, var_values)
            x_pos = np.arange(barpos,len(matrices) + barpos)
            barpos += increment
            plt.bar(x_pos,var_values[algo],label=f"{algoname} ", width = width, hatch = "//", edgecolor = "black")
        if len(var_values) == 0: 
            plt.close()    
            continue
        
        dense_amp = np.nanmean(var_values[5]/var_values[2])
        dense_amp_orig = np.nanmean(var_values[5]/res_df["density"].values)
        print("AVG AMPLIFICATION", dense_amp)
        print("AVG AMPLIFICATION, original", dense_amp_orig)
        plt.text(1, 1, f"dense-amp: {dense_amp}, {dense_amp_orig}", bbox=dict(fill=False, edgecolor='red', linewidth=2))

        x_pos = np.arange(barpos,len(matrices) + barpos)
        barpos += increment
        plt.bar(x_pos, res_df["density"].values, label = "no blocking (original)", width = width, hatch = "..", color = "white", edgecolor='black')
        

        plt.legend()
        plt.yscale("log")
        plt.title(f"block height = {row_block_size}, block width = {col_block_size}")
        plt.xticks(range(len(matrices)), matrices, rotation=90)
        plt.savefig(savename + f"_{row_block_size}_{col_block_size}.png",  bbox_inches='tight', dpi = 300)
        plt.close()    
"""




def barplot(x_labels, x_ax_label, ys, y_labels, y_styles, y_ax_label, yscale = "log", savename = "test_barplot"):
    plt.figure()
    plt.xlabel(x_ax_label)
    plt.ylabel(y_ax_label)
    
    bars = len(ys)
    tot_width = 0.9
    barpos = -tot_width/2
    increment = tot_width/bars
    width = increment*0.95

    for y, label, style in zip(ys, y_labels, y_styles):
        x_pos = np.arange(barpos,len(x_labels) + barpos)
        plt.bar(x_pos,y,label=label, width = width, **style)
        barpos += increment

    plt.legend()
    plt.yscale(yscale)
    plt.grid("both")
    plt.title(f"block height = {row_block_size}, block width = {col_block_size}")
    plt.xticks(range(len(x_labels)), x_labels, rotation=45)
    plt.savefig(savename + ".png",  bbox_inches='tight', dpi = 300)
    plt.close()    



#data_file = "test_all_suitsparse.csv"
#exp_name = "suitsparse_all"

data_file = "test_suitsparse_3.csv"
exp_name = "suitsparse_3"
image_folder = f"images/blocking_images/{exp_name}"
try: os.mkdir(image_folder) 
except: 1

df = pd.read_csv(data_file)
#df = get_dataframe_folder("results/suitsparse_collection_3")


print(df[df["VBR_nzblocks_count"] == 0])
df = df.loc[df.groupby(["matrix","blocking_algo","col_block_size","row_block_size"])["VBR_nzblocks_count"].idxmin()]

df["density"] = df["nonzeros"].values/(df["rows"].values * df["cols"].values)
df["block_density"] = df["nonzeros"].values/df["VBR_nzcount"].values
df["dense_amp"] = df["block_density"].values/df["density"].values
df["block_area"] = df["row_block_size"].values*df["col_block_size"].values
df["block ratio"] = df["row_block_size"].values/df["col_block_size"].values
df.sort_values(by=['density','matrix'], inplace=True)

#PREPARE HEATMAP FOR ALL BLOCK-SIZES

colormap_variable = "dense_amp"
max_var_value = 50


print("TESTING DF MERGE****************")
heatmap_df = pd.merge(df[df["blocking_algo"] == 2], df[df["blocking_algo"] == 5], how = "left", on = ["matrix","col_block_size","row_block_size"], suffixes=('_2','_5'))
heatmap_df["relative_dense_amp"] = heatmap_df["dense_amp_5"]/heatmap_df["dense_amp_2"]
print(heatmap_df.columns)
print(heatmap_df.loc[0].to_string())




x_var = "density"
y_var = "block_density"
c_var = "block ratio"
#MAKE SCATTER PLOT 
for blocking_algo in (2,5):
    for area in (64*64, 128*128, 256*256):
        try:
            fig,ax = plt.subplots()
            tmp_df = df.loc[(df["blocking_algo"] == blocking_algo) & (df["block_area"] == area)]
            ax = tmp_df.plot.scatter(
                            x= x_var,
                            y= y_var,
                            c= c_var,
                            colormap='viridis',alpha=0.5, edgecolor = "black")
            plt.ylabel(label_dict[y_var])
            plt.xscale("log")
            plt.yscale("log")
            plt.savefig(f"{image_folder}/{exp_name}_scatterplot_{x_var}_vs_{y_var}_area_{(area)**0.5}_algo_{blocking_algo}.png", bbox_inches='tight', dpi = 300)
        except:
            print("FAILED TO PRODUCE SCATTER PLOT FOR AREA", area)
        plt.close()


x_var = "density"
y_var = "block_density"
for row_block_size in (32,64,128,256,512,1024):
    for col_block_size in (32,64,128,256,512,1024):
        try:
            fig,ax = plt.subplots()
            tmp_df = df.loc[(df["blocking_algo"] == 5) & (df["row_block_size"] == row_block_size) & (df["col_block_size"] == row_block_size)]
            tmp_df.plot.scatter(
                            x= x_var,
                            y= y_var,
                            colormap='viridis',alpha=0.5, edgecolor = "black", marker = "*", ax = ax, label = "blocking, reordering")
            tmp_df = df.loc[(df["blocking_algo"] == 2) & (df["row_block_size"] == row_block_size) & (df["col_block_size"] == row_block_size)]
            tmp_df.plot.scatter(
                            x= x_var,
                            y= y_var,
                            colormap='viridis',alpha=0.5, edgecolor = "black", marker = "s", ax = ax, label = "blocking, no reordering")
            plt.ylabel(label_dict[y_var])
            plt.xscale("log")
            plt.yscale("log")
            plt.savefig(f"{image_folder}/{exp_name}_scatterplot_{x_var}_vs_{y_var}_{row_block_size}x{col_block_size}.png", bbox_inches='tight', dpi = 300)
        except:
            print(f"FAILED TO PRODUCE SCATTER PLOT FOR {row_block_size} x {col_block_size}")
        plt.close()



color_vars = ("relative_dense_amp","block_density_5","block_density_2","tau_5")
color_labels = ("Density Amplification (against natural blocking)", "Density Amplification (against unblocked matrix)", "Density Amplification (against unblocked matrix)","tau")

for colormap_variable,color_label in zip(color_vars, color_labels):
    vmin = 1
    vmax = 2 if (colormap_variable == "relative_dense_amp") else max_var_value
    table = heatmap_df.pivot_table(index="row_block_size", columns="col_block_size", values=colormap_variable, aggfunc='mean')
    table = table.sort_values(by=['row_block_size'], ascending=False)
    sns.heatmap(table,annot = True, cbar_kws={'label': color_label},vmin = vmin, vmax = vmax)
    plt.ylabel("block height")
    plt.xlabel("block width")
    plt.savefig(f"{image_folder}/{exp_name}_heatmap_{colormap_variable}.png",  bbox_inches='tight', dpi = 300)
    plt.close()
    

"""
for algo in (2,5):
    heatmap_df = df[df["blocking_algo"]==algo].pivot_table(index="row_block_size", columns="col_block_size", values=colormap_variable, aggfunc='mean')
    heatmap_df = heatmap_df.sort_values(by=['row_block_size'], ascending=False)
    sns.heatmap(heatmap_df,annot = True, cbar_kws={'label': color_label},vmin = 1, vmax = max_var_value)
    plt.ylabel("block height")
    plt.xlabel("block width")
    plt.savefig(f"{image_folder}/{exp_name}_heatmap_{colormap_variable}_algo_{algo}.png",  bbox_inches='tight', dpi = 300)
    plt.close()
"""




for row_block_size in sorted(df["row_block_size"].unique()):
    for col_block_size in sorted(df["col_block_size"].unique()):
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
        dense_amp = np.nanmean(data_lines["blocking, reordering"]/data_lines["blocking, no reordering"])
        dense_amp_orig = np.nanmean(data_lines["blocking, reordering"]/data_lines["no blocking (original)"])

        print(f"*******DENSE_AMP (vs blocked): {dense_amp}")
        print(f"*******DENSE_AMP (vs original): {dense_amp_orig}")

        savename = f"{image_folder}/{exp_name}_Density_barplot_{row_block_size}_{col_block_size}"

        styles = [bar_style_reordering, bar_style_blocking, bar_style_natural]

        barplot(
                x_labels = matrices_names,
                x_ax_label="graphs",
                ys=data_lines.values(),
                y_labels=data_lines.keys(),
                y_styles = styles, 
                y_ax_label = "Density", 
                savename = savename)
