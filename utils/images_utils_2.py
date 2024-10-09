import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean
import argparse
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

#____________PLOTTING PARAMS________________
color_dict = {
    "original": "#4c72b0",  # A calm, modern blue
    "clubs": "#f28e2b",     # A vibrant orange
    "metis-edge-cut": "#9467bd",  # A rich purple, distinct and balanced
    "metis-volume": "#76b7b2",     # A muted teal for balance
    "patoh": "#59a14f",    # A fresh, vibrant green
}


marker_dict = {
    "original": "o",
    "clubs": "x",
    "metis-edge-cut" : ".",
    "metis-volume" : ".",
    "patoh" : "-",
}

labels_dict = {
    "original": "Original",
    "clubs": "CLUB",
    "metis-edge-cut" : "GP-edge",
    "metis-volume" : "GP-volume",
    "patoh" : "HGP",
}

routine_labels = {
    "spmmcsr": "SpMM-CSR",
    "spmmbsr": "SpMM-BSR",
    "spmvcsr": "SpMV-CSR",
    "spmvbsr": "SpMV-BSR"

}

def count_matrices(df):
    return len(df["matrix"].unique())

def percent_improvement_formatter(y, pos):
        return f'{y*100 - 100:.0f}%'

def percent_formatter(y, pos):
        return f'{y*100:.0f}%'
#----------------------FUNCTIONS

def count_fumbles(dfs, matrices, parameter, fumbles_parameter, min_best):
    params = list(dfs.keys())
    worse_counts = {}
    for p in params:
        temp_method_df = dfs[p].copy()
        temp_method_df = set_allowed_matrices(temp_method_df, matrices)
        temp_method_df = find_best(temp_method_df, best_parameter=parameter, min=min_best)
        worse_counts[p]= (temp_method_df[fumbles_parameter] < .98).sum()
    return worse_counts
# _____________________________________________________
def count_best_method(dfs, matrices, parameter = "time_spmmcsr", min = True):
    all_results_df = pd.DataFrame()
    params_values = list(dfs.keys())

    for param in params_values:
        df = dfs[param].copy()
        df = set_allowed_matrices(df, matrices)
        df = find_best(df, best_parameter=parameter, min=min)
        df["param"] = param
        all_results_df = pd.concat([all_results_df,df])

    if "original" in params_values:
        all_results_df = all_results_df.sort_values(by="param", ascending=False, key=lambda x: x == "original") #ensures ties are win by original
    overall_best_results = find_best(all_results_df, best_parameter = parameter)

    best_param_counts = {}
    for param in params_values:
        best_param_counts[param] = (overall_best_results['param'] == param).sum()

    return best_param_counts

def calculate_improvement(best_dfs, methods, parameter = "time"):
    common_matrices = find_common_matrices(best_dfs,methods)
    improvements = {}
    for method in methods:
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df,common_matrices)
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        merged = pd.merge(original_df[['matrix', parameter]], method_df[['matrix', parameter]], on='matrix', suffixes=('_original', '_method'))
        merged['ratio'] = merged[parameter + '_original'] / merged[parameter + '_method']
        improvements[method] = gmean(merged['ratio'].values)  
    return improvements


def calculate_total_times(best_dfs, methods):
    sum_dict = {}
    common_matrices= find_common_matrices(best_dfs, methods)
    for method in methods:
        df = best_dfs[method].copy()
        clean_df = set_allowed_matrices(df,common_matrices)
        sum_dict[method] = sum(clean_df["time"])

    return sum_dict

def calculate_total_best_time(best_dfs, methods):
    common_matrices= find_common_matrices(best_dfs, methods)
    all_results_df = pd.DataFrame()
    for method in methods:
        df = best_dfs[method].copy()
        df = set_allowed_matrices(df,common_matrices)
        df["method"] =  method
        all_results_df = pd.concat([all_results_df,df])
    
    all_results_df = all_results_df.sort_values(by="method", ascending=False, key=lambda x: x == "original") #ensures ties are win by original
    overall_best_results = find_best(all_results_df)
    
    return sum(overall_best_results["time"]) 

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

def find_best(df, params = [], best_parameter = "time", min = True):
    df = df.reset_index(drop=True)  # Ensure unique index
    params = ["matrix",] + params

    if min:
        df[best_parameter] = df[best_parameter].fillna(float('inf'))
        return df.loc[df.groupby(params)[best_parameter].idxmin()]
    else:
        df[best_parameter] = df[best_parameter].fillna(float('-inf'))
        return df.loc[df.groupby(params)[best_parameter].idxmax()]


def set_allowed_matrices(df, allowed_matrices):
    return df[(df["matrix"]).isin(allowed_matrices)]

def find_common_matrices(dfs, methods):
    matrices_with_original = dfs["original"]['matrix'].unique()
    common_matrices=set(matrices_with_original)
    for method in methods:
        if method != "original":
            common_matrices &= set(dfs[method]['matrix'].unique())
    return common_matrices


def reordering_time_comparison(df,df_reordering_times):

   common_columns = list(set(df.columns) & set(df_reordering_times.columns))
   merged_df = df_reordering_times.merge(
                df[common_columns + ["time_spmmcsr","time_spmmbsr","time_spmmcsr_original","time_spmmbsr_original"]],
                on=list(common_columns),
                how="inner",
            ) 
      
   # Initialize bins for number_recover
   bins = [0, 10**3, 10**5, 10**6,10**7,10**20]
   bin_labels = ['0-10^3', '10^3-10^5', '10^5-10^6','10^6-10^7', ">10^7"]

   for routine in ["spmmcsr","spmmbsr"]:
       print("ROUTINE:",routine)
       merged_df[f"reorder_ratio_{routine}"] = merged_df["reordering_time"]/merged_df[f"time_{routine}_original"]
       merged_df[f"number_recover_{routine}"] = merged_df["reordering_time"]/(merged_df[f"time_{routine}_original"] - merged_df[f"time_{routine}"])
       merged_df = find_best(merged_df,best_parameter="reordering_time")
       good_df = merged_df.loc[merged_df[f"time_{routine}_original"] > merged_df[f"time_{routine}"]].copy()
       good_df.loc[:,f"number_recover_bin_{routine}"] = pd.cut(good_df[f"number_recover_{routine}"], bins=bins, labels=bin_labels, include_lowest=True)

       # Count the number of occurrences in each bin and calculate the percentage
       counts = good_df[f"number_recover_bin_{routine}"].value_counts(normalize=True) * 100
       counts = counts.reindex(bin_labels, fill_value=0)  # Ensure all bins are represented even if they are empty

       # Display the result
       print(f"Percentage of matrices for {routine} in each bin:")
       print(counts)


def best_barplot(dfs_reordering, square_matrices, rectangular_matrices, methods, fumbles = True, fumbles_parameter = "speedup_spmmcsr", min_best = True, parameter = "time_spmmcsr", ylabel = "", save_path = ""):
    fumbles_colors = ["gray"]
    square_counts = count_best_method(dfs_reordering, square_matrices, parameter)
    rect_counts = count_best_method(dfs_reordering, rectangular_matrices, parameter)
    #worse_counts = count_fumbles(df, )

    worse_counts = count_fumbles(dfs_reordering, 
                                matrices=square_matrices, 
                                parameter = parameter,
                                fumbles_parameter = fumbles_parameter, 
                                min_best = min_best)
    print(f"FUMBLES for {save_path}: {worse_counts}")
    for method in methods:
        temp_method_df = dfs_reordering[method].copy()
        temp_method_df = set_allowed_matrices(temp_method_df, square_matrices)
        temp_method_df = find_best(temp_method_df, best_parameter=parameter, min=min_best)
        worse_counts[method]= (temp_method_df[fumbles_parameter] < .98).sum()
    worse_counts = [-worse_counts[val] for val in methods]



    square_number = sum(square_counts.values())
    rect_number = sum(rect_counts.values())

    methods_names = [labels_dict[method] for method in methods]
    fig, ax = plt.subplots(figsize=(10, 6))
    bars_square = plt.bar(methods_names, list(square_counts.values()), 
                          label='Square Matrices',
                          edgecolor = "black", 
                          color=[color_dict[method] for method in methods])

    # Plot rectangular matrices on top of square matrices (second layer)
    bars_rect = plt.bar(methods_names, list(rect_counts.values()), 
                        label='Rectangular Matrices', 
                        bottom=list(square_counts.values()), 
                        hatch='//', 
                        edgecolor = "black",
                        color=[color_dict[method] for method in methods])

    
    if fumbles: plt.bar(methods_names, worse_counts,
                label='Fumbles Count (Worse Than Original)', 
                color=fumbles_colors, 
                edgecolor="black")


    # Shade the bottom part of the graph (from ymin to 0) in light red


    plt.ylabel(ylabel)


    #global legend
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=len(methods), frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed

    legend_elements = [Patch(facecolor='white', edgecolor='grey', label=f'Square Matrices ({square_number})'),
                    Patch(facecolor='white', edgecolor='grey', label=f'Rectangular Matrices ({rect_number})', hatch = "//"),
                    Patch(facecolor=fumbles_colors[0], edgecolor='grey', label=f'Fumbles')]

    # Add custom legend
    plt.legend(handles=legend_elements)
    
    
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12, length=6, width=2, direction='inout')
    plt.tick_params(axis='both', which='minor', length=4, width=1, direction='inout')

    plt.draw()  # Ensure the plot is fully drawn to get correct axis limits

    current_ymin, current_ymax = ax.get_ylim()
    ax.axhspan(ymin=current_ymin, ymax=0, facecolor='#FFCCCC', alpha=0.2)


    plt.tight_layout()

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def best_barplot_parameter(dfs_reordering,square_matrices, rectangular_matrices, method, plot_params_values = [], plot_parameter = "mask",improvement_parameter = "None", ratio_parameter = "None", xlabel = "", ylabel = "", save_path = "", min_best = True):

    df = dfs_reordering[method].copy()
    colors = ['#FFD699', '#FFB84D', '#FF8C00', '#CC7000']
    fumbles_colors = ["gray"]
    #separate df depending on parameter value
    if not plot_params_values: plot_params_values = df[plot_parameter].unique()
    plot_params_labels = [f"{p}" for p in plot_params_values]
    dfs_parameter = {}
    worse_counts = {}
    for plot_parameter_value in plot_params_values:
        dfs_parameter[plot_parameter_value] = df.loc[df[plot_parameter]==plot_parameter_value]
        dfs_parameter[plot_parameter_value] = find_best(dfs_parameter[plot_parameter_value], best_parameter=improvement_parameter, min=min_best)
        worse_counts[plot_parameter_value] = (dfs_parameter[plot_parameter_value][ratio_parameter] < .98).sum()
        print(worse_counts[plot_parameter_value])
    square_counts = count_best_method(dfs_parameter, square_matrices,improvement_parameter)
    rect_counts = count_best_method(dfs_parameter, rectangular_matrices,improvement_parameter)

    print(square_counts, rect_counts)
    square_number = sum(square_counts.values())
    rect_number = sum(rect_counts.values())

    fig, ax = plt.subplots(figsize=(10, 6))
    bars_square = plt.bar(plot_params_labels, list(square_counts.values()), 
                          label='Square Matrices',
                          edgecolor = "black",
                          color = colors)

    # Plot rectangular matrices on top of square matrices (second layer)
    bars_rect = plt.bar(plot_params_labels, list(rect_counts.values()), 
                        label='Rectangular Matrices', 
                        bottom=list(square_counts.values()), 
                        hatch='//',
                        color = colors, 
                        edgecolor = "black")

    plt.bar(plot_params_labels, [-worse_counts[val] for val in plot_params_values],
                label='Fumbles Count (Worse Than Original)', 
                color=fumbles_colors, edgecolor="black")



    plt.ylabel(ylabel)
    plt.xlabel(ylabel)


    legend_elements = [Patch(facecolor='white', edgecolor='grey', label=f'Square Matrices ({square_number})'),
                    Patch(facecolor='white', edgecolor='grey', label=f'Rectangular Matrices ({rect_number})', hatch = "//"),
                    Patch(facecolor=fumbles_colors[0], edgecolor='grey', label=f'Fumbles')]

    # Add custom legend
    plt.legend(handles=legend_elements)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12, length=6, width=2, direction='inout')
    plt.tick_params(axis='both', which='minor', length=4, width=1, direction='inout')


    plt.draw()  # Ensure the plot is fully drawn to get correct axis limits

    current_ymin, current_ymax = ax.get_ylim()
    ax.axhspan(ymin=current_ymin, ymax=0, facecolor='#FFCCCC', alpha=0.2)


    plt.tight_layout()

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def make_improvements_barplot_and_distribution(dfs_reordering, methods, matrices, set_negative_1=False, set_missing_1=False, min_best=False, parameter="None", ylabel="Speedup", save_path="test.png", verbose=False):
    """
    Creates a barplot with error bars and a sideways histogram overlayed on the right.
    """

    methods = methods.copy()
    if "original" in methods:
        methods.remove("original")

    improvements = {}
    low_errors = []
    up_errors = []
    hist_data = {}  # To store data for histogram
    
    all_values = []  # Collect all data to determine common bins
    for method in methods:
        df = dfs_reordering[method].copy()
        df = set_allowed_matrices(df, matrices)
        df = find_best(df, best_parameter=parameter, min=min_best)

        if not set_missing_1 and df[parameter].isna().any():
            print(f"ERROR: MISSING VALUES for {parameter}")

        if set_missing_1:
            df.loc[df[parameter].isna()] = 1

        if set_negative_1:
            df.loc[df[parameter] < 1] = 1

        # Collect data for histogram
        hist_data[method] = df[parameter].values
        all_values.extend(df[parameter].values)  # Add data to determine common bins

        # Calculate percentiles for the bar plot
        percentiles = np.percentile(df[parameter].values, [25, 50, 75])
        improvements[method] = percentiles[1] 
        if verbose: 
            print(method, percentiles, improvements[method])
        low_errors.append(improvements[method] - percentiles[0])
        up_errors.append(percentiles[2] - improvements[method])

    # Calculate appropriate y-limits
    max_y = (max(up_errors) + max(improvements.values())) * 1.01
    min_y = min(1, (min(improvements.values()) - min(low_errors)) * 0.9)
    
    # Determine common bin edges for all histograms
    #common_bins = [np.histogram_bin_edges(all_values, bins=200)]
    common_bins = np.linspace(min_y,max_y,20)

    # Create figure and gridspec to manage layout
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)

    # Plot the bar chart
    ax_point = fig.add_subplot(gs[0])
    for i, method in enumerate(methods):
        method_name = labels_dict[method]
        ax_point.errorbar(method_name, improvements[method], yerr=[[low_errors[i]], [up_errors[i]]],
                          color=color_dict[method], ecolor=color_dict[method], 
                          elinewidth=2, capsize=5, marker='o', markersize=10,label=method)


    ax_point.set_xlabel("Reordering technique")
    ax_point.set_ylabel(ylabel)  # Set ylabel to "Speedup"
    ax_point.set_ylim(min_y, max_y)
 
    # Format y-axis as percentage
    def speedup_formatter(y, pos):
        return f'{y*100 - 100:.0f}%'
    
    ax_point.yaxis.set_major_formatter(FuncFormatter(speedup_formatter))
    
    # Ensure that the y-ticks are visible for the bar plot
    ax_point.tick_params(axis='y', which='both', labelleft=True)

    # Plot the sideways histogram without sharey
    ax_hist = fig.add_subplot(gs[1])

    offsets = np.linspace(-0.05, 0.05, len(methods)) * (common_bins[1] - common_bins[0])

    for i, method in enumerate(methods):
        # Slightly shift each histogram by a small offset to avoid overlap
        ax_hist.hist(hist_data[method], bins=common_bins + offsets[i], orientation='horizontal', 
                     alpha=0.6, color=color_dict[method], histtype='step', label=method, linewidth=2)

    ax_hist.set_xlabel("Frequency")
    ax_hist.set_yticklabels([])  # Hide the y-ticks on the histogram
    ax_hist.set_xlim(0, None)  # Auto-adjust x-limits for the histogram

    # Sync the y-limits of the histogram with the bar chart
    ax_hist.set_ylim(ax_point.get_ylim())

    # Adjust layout and save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def make_improvements_barplot(dfs_reordering, methods, matrices, allow_missing=False, min_best=False, parameter="None", ylabel="Speedup", save_path="test.png", verbose=False):
    """
    Creates a barplot with error bars and a sideways histogram overlayed on the right.
    """
    methods = methods.copy()
    if "original" in methods:
        methods.remove("original")

    print("MAKING IMP PLOT")
    improvements = {}
    low_errors = []
    up_errors = []
    hist_data = {}  # To store data for histogram
    
    all_values = []  # Collect all data to determine common bins
    for method in methods:
        df = dfs_reordering[method].copy()
        df = set_allowed_matrices(df, matrices)
        df = find_best(df, best_parameter=parameter, min=min_best)

        if not allow_missing and df[parameter].isna().any():
            print(f"ERROR: MISSING VALUES for {parameter}")

        if allow_missing:
            df.loc[df[parameter].isna(), parameter] = float("-inf")

        df.loc[df[parameter] == float("-inf") , parameter] = 0
        df.loc[df[parameter] == float("inf") , parameter] = 10

        # Collect data for histogram
        hist_data[method] = df[parameter].values
        all_values.extend(df[parameter].values)  # Add data to determine common bins

        # Calculate percentiles for the bar plot
        percentiles = np.percentile(df[parameter].values, [25, 50, 75])
        print(method, percentiles)
        improvements[method] = percentiles[1] 
        if verbose: 
            print(method, percentiles, improvements[method])
        low_errors.append(improvements[method] - percentiles[0])
        up_errors.append(percentiles[2] - improvements[method])

    # Calculate appropriate y-limits
    max_y = (max(up_errors) + max(improvements.values())) * 1.01
    min_y = min(1, (min(improvements.values()) - min(low_errors)) * 0.9)
    
    fig = plt.figure(figsize=(10, 6))
    
    # Plot the bar chart
    for i, method in enumerate(methods):
        method_name = labels_dict[method]
        plt.errorbar(method_name, improvements[method], yerr=[[low_errors[i]], [up_errors[i]]],
                          color=color_dict[method], ecolor=color_dict[method], 
                          elinewidth=2, capsize=5, marker='o', markersize=10,label=labels_dict[method])

    #ax_point.legend(loc='best', frameon=False)
    plt.gca().set_ylabel(ylabel)  # Set ylabel to "Speedup"
    plt.gca().set_ylim(min_y, max_y)
    plt.gca().grid(True, linestyle='--', alpha=0.3)
    plt.gca().tick_params(axis='both', which='major', length=5, direction='in')
    
    plt.gca().yaxis.set_major_formatter(FuncFormatter(percent_improvement_formatter))
    
    # Ensure that the y-ticks are visible for the bar plot
    plt.gca().tick_params(axis='y', which='both', labelleft=True)
    
    #global legend
    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=len(methods), frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed


    # Adjust layout and save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def make_improvements_barplot_and_distribution_2(dfs_reordering, methods, matrices, allow_missing=False, min_best=False, parameter="None", ylabel="Speedup", save_path="test.png", verbose=False):
    """
    Creates a barplot with error bars and a sideways histogram overlayed on the right.
    """
    methods = methods.copy()
    if "original" in methods:
        methods.remove("original")

    print("MAKING IMP PLOT")
    improvements = {}
    low_errors = []
    up_errors = []
    hist_data = {}  # To store data for histogram
    
    all_values = []  # Collect all data to determine common bins
    for method in methods:
        df = dfs_reordering[method].copy()
        df = set_allowed_matrices(df, matrices)
        df = find_best(df, best_parameter=parameter, min=min_best)

        if not allow_missing and df[parameter].isna().any():
            print(f"ERROR: MISSING VALUES for {parameter}")

        if allow_missing:
            df.loc[df[parameter].isna(), parameter] = float("-inf")

        df.loc[df[parameter] == float("-inf") , parameter] = 0
        df.loc[df[parameter] == float("inf") , parameter] = 10

        # Collect data for histogram
        hist_data[method] = df[parameter].values
        all_values.extend(df[parameter].values)  # Add data to determine common bins

        # Calculate percentiles for the bar plot
        percentiles = np.percentile(df[parameter].values, [25, 50, 75])
        print(method, percentiles)
        improvements[method] = percentiles[1] 
        if verbose: 
            print(method, percentiles, improvements[method])
        low_errors.append(improvements[method] - percentiles[0])
        up_errors.append(percentiles[2] - improvements[method])

    # Calculate appropriate y-limits
    max_y = (max(up_errors) + max(improvements.values())) * 1.01
    min_y = min(1, (min(improvements.values()) - min(low_errors)) * 0.9)
    
    # Determine common bin edges for all histograms
    #common_bins = [np.histogram_bin_edges(all_values, bins=200)]
    common_bins = np.linspace(min_y,max_y,20)

    # Create figure and gridspec to manage layout
    fig = plt.figure(figsize=(10, 6))
    figures = len(methods) + 1
    gs = fig.add_gridspec(1, figures, width_ratios=[2,] + [1]*len(methods), wspace=0.05)


    # Plot the bar chart
    ax_point = fig.add_subplot(gs[0])
    for i, method in enumerate(methods):
        method_name = labels_dict[method]
        ax_point.errorbar(method_name, improvements[method], yerr=[[low_errors[i]], [up_errors[i]]],
                          color=color_dict[method], ecolor=color_dict[method], 
                          elinewidth=2, capsize=5, marker='o', markersize=10,label=labels_dict[method])
    ax_point.set_xticklabels([])
    ax_point.set_xlabel("Reordering Technique")
    #ax_point.legend(loc='best', frameon=False)
    ax_point.set_ylabel(ylabel)  # Set ylabel to "Speedup"
    ax_point.set_ylim(min_y, max_y)
    ax_point.grid(True, linestyle='--', alpha=0.3)
    ax_point.tick_params(axis='both', which='major', length=5, direction='in')
    
    ax_point.yaxis.set_major_formatter(FuncFormatter(percent_improvement_formatter))
    
    # Ensure that the y-ticks are visible for the bar plot
    ax_point.tick_params(axis='y', which='both', labelleft=True)
    
    # Plot the sideways histogram without sharey

    def hist_formatter(y, pos):
        return f'{y*100/len(matrices):.0f}%'

    hist_axes = [None for i in range(len(methods))]
    for i, method in enumerate(methods):
        hist_axes[i] = fig.add_subplot(gs[1 + i])
        hist_axes[i].hist(hist_data[method], bins=common_bins, orientation='horizontal', 
                     alpha=0.9, color=color_dict[method], histtype='step', label=method, linewidth=2)
        hist_axes[i].set_xlim(0, 65)

        # Sync the y-limits of the histogram with the bar chart
        hist_axes[i].set_ylim(ax_point.get_ylim())
        #hist_axes[i].xaxis.set_major_formatter(FuncFormatter(hist_formatter))
        hist_axes[i].set_yticklabels([])  # Hide the y-ticks on the histogram
        hist_axes[i].xaxis.set_major_locator(MaxNLocator(nbins=4))
        hist_axes[i].grid(True, linestyle='--', alpha=0.3)
        hist_axes[i].tick_params(axis='x', labelsize = 9)  


    #make the x-label for the histograms
    left = hist_axes[0].get_position().x0
    right = hist_axes[-1].get_position().x1
    bottom = hist_axes[0].get_position().y0*1.08  # Assuming all histograms are aligned vertically
    x_center = (left + right) / 2
    y_position = bottom - 0.06  # Adjust as needed to position the label below the histograms
    fig.text(x_center, y_position, "Frequency", ha='center', va='top')

    current_x_coord, current_y_coord = ax_point.xaxis.get_label().get_position()
    ax_point.xaxis.set_label_coords(current_x_coord, current_y_coord - 0.055 )  # Adjust y-coordinate only

    #global legend
    handles, labels = ax_point.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.07),
           ncol=len(methods), frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed



    # Adjust layout and save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_speedup_vs_matrix_param(best_dfs, methods, matrices = None, param="density", title='Speedup vs matrix density', save_path="test_hist.png"):



    common_matrices = find_common_matrices(best_dfs, methods)
    if matrices != None: 
        common_matrices &= matrices

    plt.figure(figsize=(10, 6))
    
    for method in methods:
        original_df = best_dfs["original"].copy()
        original_df = set_allowed_matrices(original_df, common_matrices)
        
        # Include the 'param' column in the merge operation
        method_df = best_dfs[method].copy()
        method_df = set_allowed_matrices(method_df, common_matrices)
        
        # Merge including the 'param' column
        merged = pd.merge(original_df[['matrix', param, 'time']], 
                          method_df[['matrix', 'time']], 
                          on=['matrix'], 
                          suffixes=('_original', '_method'))
        
        merged = merged.sort_values(by=param)
        
        # Calculate the speedup ratio
        merged['ratio'] = merged['time_method'] / merged['time_original']
        speedups = 1 / merged['ratio']
        
        # Ensure speedups are capped at 1
        speedups[speedups < 1] = 1
        
        # Extract the values of the parameter to plot against
        param_values = merged[param]
        
        # Plot the data
        plt.plot(param_values, speedups, color = color_dict[method], marker='o', linestyle='',markersize= 3, label=f'{method}')
    
    plt.xscale('log')
    plt.ylim(0.9, 1.5)
    plt.xlabel(param.capitalize())
    plt.ylabel('Speedup')
    #plt.title(title)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_improvement_by_matrix_old(dfs_reordering, order_by, methods, parameter = "None", matrices = None, min_best = False, title = "", yFormatter = percent_formatter,ylim = [0,5], xlabel = "Default x Label", ylabel = "Default Y Label", save_path="test_speedup_by_matrix.png"):
    """
    Plots the speedup for each matrix, sorted by speedups, for each method.

    Parameters:
    - best_dfs: Dictionary of DataFrames for each method containing matrix times.
    - methods: List of methods to compare.
    """

    new_df = {}

    for method in methods:
        matrices &= set(dfs_reordering[method]["matrix"].unique())

    for method in methods:
        new_df[method] = dfs_reordering[method].copy()
        new_df[method] = set_allowed_matrices(new_df[method], matrices)
        new_df[method] = find_best(new_df[method], best_parameter=parameter, min=min_best)

    # First, calculate the speedup for the 'clubs' method and store the ordering
    if not order_by in methods: 
        print("WARNING: plot_speedup_by_matrix: ORDER_BY must appear also in METHODS")
        exit(1)
        return

    plt.figure(figsize=(10, 6))

    new_df[order_by] = new_df[order_by].sort_values(by=parameter, ascending = not min_best)
    ordered_matrices = new_df[order_by]['matrix'].values
    plt.plot(range(len(ordered_matrices)), new_df[order_by][parameter], 
                  color=color_dict[order_by], 
                  marker='o', 
                  linestyle='', 
                  markersize=3, 
                  label=f'{labels_dict[method]}')
    
    for method in methods:
        if method == order_by: continue
        new_df[method]['matrix'] = pd.Categorical(new_df[method]['matrix'], categories=ordered_matrices, ordered=True)
        plt.plot(range(len(ordered_matrices)), new_df[method][parameter], 
                 color=color_dict[method], 
                 marker='o', 
                 linestyle='', 
                 markersize=3, 
                 label=f'{labels_dict[method]}')


    #plt.yscale("log")

    plt.gca().yaxis.set_major_formatter(FuncFormatter(yFormatter))

    
    plt.axhline(0, color=color_dict["original"])
    plt.ylim(ylim[0],ylim[1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.title(title)


    #top legend
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=len(methods) + 1, frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed

    plt.gca().grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_improvement_by_matrix(dfs_reordering, order_by, methods, parameter="None", matrices=None, min_best=False, title="", y_scale = "linear", original_line_y = 0, x_scale = "linear", yFormatter=lambda y, pos: y, ylim=[0, 5], xlabel="Default x Label", ylabel="Default Y Label", save_path="test_speedup_by_matrix.png"):
    """
    Plots the speedup for each matrix, sorted by speedups, for each method.

    Parameters:
    - dfs_reordering: Dictionary of DataFrames for each method containing matrix times.
    - order_by: Method used to order matrices.
    - methods: List of methods to compare.
    - parameter: Parameter to plot.
    - matrices: Set of matrices to include.
    - min_best: Boolean indicating if lower values are better.
    - yFormatter: Function to format y-axis labels.
    - ylim: List containing y-axis limits [min, max].
    - xlabel: Label for the x-axis.
    - ylabel: Label for the y-axis.
    - save_path: Path to save the plot image.
    """

    edgecolor = (0, 0, 0, 0.3)



    new_df = {}

    for method in methods:
        matrices &= set(dfs_reordering[method]["matrix"].unique())

    for method in methods:
        new_df[method] = dfs_reordering[method].copy()
        new_df[method] = set_allowed_matrices(new_df[method], matrices)
        new_df[method] = find_best(new_df[method], best_parameter=parameter, min=min_best)

    if order_by not in methods:
        print("WARNING: plot_speedup_by_matrix: ORDER_BY must appear also in METHODS")
        exit(1)
        return

    plt.figure(figsize=(10, 6))

    # Sorting and ordering matrices
    new_df[order_by] = new_df[order_by].sort_values(by=parameter, ascending=not min_best)
    ordered_matrices = new_df[order_by]['matrix'].values
    x_values = np.arange(len(ordered_matrices))

    # Setting up delta for y-axis limits
    delta = 0.05 * (ylim[1] - ylim[0])
    plt.ylim(ylim[0] - delta, ylim[1] + delta)

    # Custom y-axis formatter to handle infinite values
    def custom_yFormatter(y, pos):
        if y < ylim[0]:
            return ""
        elif y < ylim[0]*1.1:
            return 'Failed'
        elif y > ylim[1]:
            return ""
        else:
            return yFormatter(y, pos)

    plt.gca().yaxis.set_major_formatter(FuncFormatter(custom_yFormatter))

    # Function to plot data for a method
    def plot_method_data(method, label_method=True):
        y_values = new_df[method][parameter].copy()
        neg_inf_mask = np.isneginf(y_values)
        pos_inf_mask = np.isposinf(y_values)
        finite_mask = np.isfinite(y_values)

        y_values_adjusted = y_values.copy()
        y_values_adjusted[neg_inf_mask] = ylim[0]
        y_values_adjusted[pos_inf_mask] = ylim[1]

        # Plot finite values
        plt.scatter(x_values[finite_mask], y_values_adjusted[finite_mask],
                 color=color_dict[method],
                 marker='o',
                 alpha=0.7,
                 linestyle='',
                 edgecolor=edgecolor,
                 linewidths=0.1,
                 s=35,
                 label=f'{labels_dict[method]}' if label_method else "_nolegend_")

        # Plot -inf values
        plt.scatter(x_values[neg_inf_mask], y_values_adjusted[neg_inf_mask],
                 color=color_dict[method],
                 marker='x',
                 s=35,
                 label='_nolegend_')

        # Plot +inf values
        if False:
            plt.scatter(x_values[pos_inf_mask], y_values_adjusted[pos_inf_mask],
                    color=color_dict[method],
                    marker='^',
                    markersize=5,
                    label='_nolegend_')

    # Plot data for each method
    plot_method_data(order_by)
    for method in methods:
        if method == order_by:
            continue
        new_df[method]['matrix'] = pd.Categorical(new_df[method]['matrix'], categories=ordered_matrices, ordered=True)
        plot_method_data(method)

    plt.axhline(original_line_y, color=color_dict.get("original", 'black'))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Adjust y-axis ticks to include ylim[0] and ylim[1]
    ax = plt.gca()
    yticks = ax.get_yticks()
    if not any(yticks <= ylim[0]):
        yticks = np.append(yticks, ylim[0])
    if not any(yticks >= ylim[1]):
        yticks = np.append(yticks, ylim[1])
    yticks.sort()
    ax.set_yticks(yticks)

    # Top legend
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.1),
               ncol=len(methods) + 1, frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed

    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12, length=6, width=2, direction='inout')
    plt.tick_params(axis='both', which='minor', length=4, width=1, direction='inout')
    #plt.yscale(y_scale)
    #plt.xscale(x_scale)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def speedup_vs_nnz_ratio(dfs_reorderings, method, x_parameter, y_parameter, matrices, xlim = [0,2], ylim=[0, 2], xscale = "linear", yscale="linear",xlabel="Default x Label", ylabel="Default Y Label", save_path="test_ratio_study.png"):
    edgecolor = (0, 0, 0, 0.3)
    
    
    df = dfs_reorderings[method].copy()
    df = set_allowed_matrices(df, matrices)
    df = find_best(df, best_parameter=x_parameter, min=False)
    df = df.sort_values(by=x_parameter)

    fig = plt.figure(figsize=(10, 6))
    plt.scatter(df[x_parameter].values, df[y_parameter].values, marker = "o", color = color_dict[method], s=50, linewidths=0.1, alpha=0.7, edgecolor=edgecolor)
    plt.gca().grid(True, linestyle='--', alpha=0.3)
    plt.yscale(xscale)
    plt.xscale(yscale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)  
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12, length=6, width=2, direction='inout')
    plt.tick_params(axis='both', which='minor', length=4, width=1, direction='inout')



    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    return

def plot_improvement_by_parameter(dfs_reorderings, plot_parameter, method = "clubs", improvement_parameter="None", allow_missing = True, matrices=None, min_best=False, title="", ylim=[0, 5], xlabel="Default x Label", ylabel="Default Y Label", save_path="test_parameter_study.png"):

    df = dfs_reorderings[method].copy()
    df = set_allowed_matrices(df,matrices)
    
    #separate df depending on parameter value
    plot_params_values = df[plot_parameter].unique()
    dfs_parameter = {}
    for plot_parameter_value in plot_params_values:
        dfs_parameter[plot_parameter_value] = df.loc[df[plot_parameter]==plot_parameter_value]
        dfs_parameter[plot_parameter_value] = find_best(dfs_parameter[plot_parameter_value], best_parameter=improvement_parameter, min=min_best)
        
    #calculate common matrices
    common_matrices = matrices
    for plot_parameter_value in plot_params_values:
        df = dfs_parameter[plot_parameter_value]
        common_matrices &= set(df[~df[improvement_parameter].isna()]["matrix"].unique())

    for plot_parameter_value in plot_params_values:
        dfs_parameter[plot_parameter_value] = set_allowed_matrices(dfs_parameter[plot_parameter_value],common_matrices)


    improvements = {}
    low_errors = []
    up_errors = []
    hist_data = {}  # To store data for histogram  
    all_values = []  # Collect all data to determine common bins
    for plot_parameter_value in plot_params_values:
        df = dfs_parameter[plot_parameter_value]
        if not allow_missing and df[improvement_parameter].isna().any():
            print(f"ERROR: MISSING VALUES for {improvement_parameter}")

        df.loc[df[improvement_parameter] == float("-inf") , improvement_parameter] = 0
        df.loc[df[improvement_parameter] == float("inf") , improvement_parameter] = 10

        # Collect data for histogram
        hist_data[plot_parameter_value] = df[improvement_parameter].values
        all_values.extend(df[improvement_parameter].values)  # Add data to determine common bins

        # Calculate percentiles for the bar plot
        percentiles = np.percentile(df[improvement_parameter].values, [25, 50, 75])
        improvements[plot_parameter_value] = percentiles[1] 
        low_errors.append(improvements[plot_parameter_value] - percentiles[0])
        up_errors.append(percentiles[2] - improvements[plot_parameter_value])
        print(plot_parameter,plot_parameter_value, percentiles, improvements[plot_parameter_value])


    # Calculate appropriate y-limits
    max_y = (max(up_errors) + max(improvements.values())) * 1.01
    min_y = min(1, (min(improvements.values()) - min(low_errors)) * 0.9)
    
    # Determine common bin edges for all histograms
    #common_bins = [np.histogram_bin_edges(all_values, bins=200)]
    common_bins = np.linspace(min_y,max_y,20)

    # Create figure and gridspec to manage layout
    fig = plt.figure(figsize=(10, 6))

    # Plot the bar chart
    ax_point = fig.add_subplot(111)
    for i, plot_parameter_value in enumerate(plot_params_values):
        ax_point.errorbar(f"{plot_parameter_value}", improvements[plot_parameter_value], yerr=[[low_errors[i]], [up_errors[i]]], 
                          elinewidth=2, capsize=5, marker='o', markersize=10,label=plot_parameter_value)

    #ax_point.legend(loc='best', frameon=False)
    ax_point.set_ylabel(ylabel)  # Set ylabel to "Speedup"
    ax_point.set_ylim(min_y, max_y)
    ax_point.grid(True, linestyle='--', alpha=0.3)
    ax_point.tick_params(axis='both', which='major', length=5, direction='in')
    
    ax_point.yaxis.set_major_formatter(FuncFormatter(percent_improvement_formatter))
    
    # Ensure that the y-ticks are visible for the bar plot
    ax_point.tick_params(axis='y', which='both', labelleft=True)
    
    #global legend
    handles, labels = ax_point.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=len(plot_params_values), frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed

    # Adjust layout and save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_improvement_by_parameter_and_distribution(dfs_reorderings, plot_parameter, method = "clubs", improvement_parameter="None", allow_missing = True, matrices=None, min_best=False, title="", ylim=[0, 5], xlabel="Default x Label", ylabel="Default Y Label", save_path="test_parameter_study.png"):

    df = dfs_reorderings[method].copy()
    df = set_allowed_matrices(df,matrices)
    
    #separate df depending on parameter value
    plot_params_values = df[plot_parameter].unique()
    dfs_parameter = {}
    for plot_parameter_value in plot_params_values:
        dfs_parameter[plot_parameter_value] = df.loc[df[plot_parameter]==plot_parameter_value]
        dfs_parameter[plot_parameter_value] = find_best(dfs_parameter[plot_parameter_value], best_parameter=improvement_parameter, min=min_best)
        


    #calculate common matrices
    common_matrices = matrices
    for plot_parameter_value in plot_params_values:
        df = dfs_parameter[plot_parameter_value]
        common_matrices &= set(df[~df[improvement_parameter].isna()]["matrix"].unique())

    for plot_parameter_value in plot_params_values:
        dfs_parameter[plot_parameter_value] = set_allowed_matrices(dfs_parameter[plot_parameter_value],common_matrices)


    improvements = {}
    low_errors = []
    up_errors = []
    hist_data = {}  # To store data for histogram  
    all_values = []  # Collect all data to determine common bins
    for plot_parameter_value in plot_params_values:
        df = dfs_parameter[plot_parameter_value]
        if not allow_missing and df[improvement_parameter].isna().any():
            print(f"ERROR: MISSING VALUES for {improvement_parameter}")

        df.loc[df[improvement_parameter] == float("-inf") , improvement_parameter] = 0
        df.loc[df[improvement_parameter] == float("inf") , improvement_parameter] = 10

        # Collect data for histogram
        hist_data[plot_parameter_value] = df[improvement_parameter].values
        all_values.extend(df[improvement_parameter].values)  # Add data to determine common bins

        # Calculate percentiles for the bar plot
        percentiles = np.percentile(df[improvement_parameter].values, [25, 50, 75])
        improvements[plot_parameter_value] = percentiles[1] 
        low_errors.append(improvements[plot_parameter_value] - percentiles[0])
        up_errors.append(percentiles[2] - improvements[plot_parameter_value])
        print(plot_parameter,plot_parameter_value, percentiles, improvements[plot_parameter_value])


    # Calculate appropriate y-limits
    max_y = (max(up_errors) + max(improvements.values())) * 1.01
    min_y = min(1, (min(improvements.values()) - min(low_errors)) * 0.9)
    
    # Determine common bin edges for all histograms
    #common_bins = [np.histogram_bin_edges(all_values, bins=200)]
    common_bins = np.linspace(min_y,max_y,20)

    # Create figure and gridspec to manage layout
    fig = plt.figure(figsize=(10, 6))
    figures = len(plot_params_values) + 1
    gs = fig.add_gridspec(1, figures, width_ratios=[2,] + [1]*len(plot_params_values), wspace=0.05)


    # Plot the bar chart
    ax_point = fig.add_subplot(gs[0])
    for i, plot_parameter_value in enumerate(plot_params_values):
        ax_point.errorbar(f"{plot_parameter_value}", improvements[plot_parameter_value], yerr=[[low_errors[i]], [up_errors[i]]], 
                          elinewidth=2, capsize=5, marker='o', markersize=10,label=plot_parameter_value)

    #ax_point.legend(loc='best', frameon=False)
    ax_point.set_ylabel(ylabel)  # Set ylabel to "Speedup"
    ax_point.set_ylim(min_y, max_y)
    ax_point.grid(True, linestyle='--', alpha=0.3)
    ax_point.tick_params(axis='both', which='major', length=5, direction='in')
    
    ax_point.yaxis.set_major_formatter(FuncFormatter(percent_improvement_formatter))
    
    # Ensure that the y-ticks are visible for the bar plot
    ax_point.tick_params(axis='y', which='both', labelleft=True)
    
    # Plot the sideways histogram without sharey

    def hist_formatter(y, pos):
        return f'{y*100/len(matrices):.0f}%'

    hist_axes = [None for i in range(len(plot_params_values))]
    for i, plot_parameter_value in enumerate(plot_params_values):
        hist_axes[i] = fig.add_subplot(gs[1 + i])
        hist_axes[i].hist(hist_data[plot_parameter_value], bins=common_bins, orientation='horizontal', 
                     alpha=0.9, histtype='step', label=plot_parameter_value, linewidth=2)
        #hist_axes[i].set_xlabel("Frequency")        
        hist_axes[i].set_xlim(0, 65)  # Auto-adjust x-limits for the histogram

        # Sync the y-limits of the histogram with the bar chart
        hist_axes[i].set_ylim(ax_point.get_ylim())
        #hist_axes[i].xaxis.set_major_formatter(FuncFormatter(hist_formatter))
        hist_axes[i].set_yticklabels([])  # Hide the y-ticks on the histogram
        hist_axes[i].xaxis.set_major_locator(MaxNLocator(nbins=4))
        hist_axes[i].grid(True, linestyle='--', alpha=0.3)



    #make the x-label for the histograms
    left = hist_axes[0].get_position().x0
    right = hist_axes[-1].get_position().x1
    bottom = hist_axes[0].get_position().y0*1.08  # Assuming all histograms are aligned vertically
    x_center = (left + right) / 2
    y_position = bottom - 0.05  # Adjust as needed to position the label below the histograms
    fig.text(x_center, y_position, "Frequency", ha='center', va='top')


    #global legend
    handles, labels = ax_point.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),
           ncol=len(plot_params_values), frameon=False)
    plt.subplots_adjust(top=1)  # Adjust the top margin as needed



    # Adjust layout and save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()



def plot_club_performance_by_param(dfs, method, performance_param = "time", param = "mask", title='Performance Variation under Mask', save_name=""):
    """
    Plots the performance variation for a specific method under different mask values,
    with one curve for each tau value, keeping centroid fixed.
    """
    # Filter the DataFrame for the selected method and fixed centroid value
    
    method_df = dfs[method].copy()
    original_df = dfs["original"]
    method_df = find_best(method_df, [param,], best_parameter= performance_param)
    
    plt.figure(figsize=(10, 6))

    merged = pd.merge(original_df[['matrix', performance_param]], 
                        method_df[['matrix', param,performance_param]], 
                        on=['matrix'], 
                        suffixes=('_original', '_method'))
    
    merged = merged.sort_values(by=param)

    merged['speedup'] = np.maximum(1, merged[performance_param + '_original'] / merged[performance_param + '_method'])

    # Ensure speedups are capped at 1

    plot_data = merged[[param,"speedup"]]

    #adds best plot
    best_merged = merged.loc[merged.groupby(by="matrix")["speedup"].idxmax()]
    speeds = best_merged["speedup"].values
    best_data = pd.DataFrame({param: ["Best"] * len(best_merged), "speedup": best_merged["speedup"].values})
    #---------------

    plot_data = pd.concat([plot_data, best_data], ignore_index=True)
    plot_data['speedup_percent'] = plot_data['speedup'] * 100 - 100

    sns.violinplot(x=param, y='speedup_percent', data=plot_data, inner='quartile', color=color_dict["clubs"])
# Plot the data with asymmetric error bars
    #plt.errorbar(param_values, speedups, yerr=[lower_errors, upper_errors], color=color_dict["clubs"], marker='o', linestyle='--', markersize=3, capsize=5)
    
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}%'))
    plt.ylim(-10,50)
    plt.xlabel(param)
    plt.ylabel('Geometric mean of the speedup (%)')
    #plt.title(title)
    plt.tight_layout()
    plt.savefig(save_name)
    plt.close()


    return 
    #lineplot
    speedups = []
    lower_errors = []
    upper_errors = []

    param_values = merged[param].values
    for param_value in param_values:
        speeds = merged[merged[param] == param_value]["speedup"].values
        speedup = gmean(speeds)
        speedups.append(speedup)
        
        # Calculate percentiles for asymmetric error bars
        lower_error = speedup - np.percentile(speeds, 25)
        upper_error = np.percentile(speeds, 75) - speedup
        # Ensure errors are non-negative
        lower_errors.append(max(0, lower_error))
        upper_errors.append(max(0, upper_error))

#-____________________ END OF FUNCTIONS_____________________________

