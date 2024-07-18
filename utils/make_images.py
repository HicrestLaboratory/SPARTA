import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define the directory containing the files
directory = 'results/results_2024'
output_plot_dir = 'results/results_2024/plots'
os.makedirs(output_plot_dir, exist_ok=True)


# Initialize empty dataframes
clubs_df = pd.DataFrame()
metis_df = pd.DataFrame()

# Column names for the dataframes

original_columns = [
    "matrix_name", "rows", "cols", "nnz", "block_size", "scramble",
    "VBR_nzcount", "VBR_nzblocks_count", "VBR_average_height", "VBR_longest_row"
]

clubs_columns = [
    "matrix_name", "rows", "cols", "nnz", "block_size", "scramble", "algo", "mask", "tau", "centroid",
    "VBR_nzcount", "VBR_nzblocks_count", "VBR_average_height", "VBR_longest_row"
]
metis_columns = [
    "matrix_name", "rows", "cols", "nnz", "block_size", "scramble", "algo", "metis_obj", "metis_part",
    "VBR_nzcount", "VBR_nzblocks_count", "VBR_average_height", "VBR_longest_row"
]

# Read and concatenate the data into dataframes
clubs_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith('clubs')]
metis_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith('metis')]
original_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith('original')]


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

clubs_df = read_and_concat(clubs_files)
metis_df = read_and_concat(metis_files)
original_df = read_and_concat(original_files)

# Function to get the best result for each matrix
def get_best_results(df, group_cols, value_col='VBR_nzblocks_count'):
    return df.loc[df.groupby(group_cols)[value_col].idxmin()]

# Function to create and save bar plots
def create_bar_plots(original_df, clubs_df, metis_df, output_dir):
    block_sizes = clubs_df['block_size'].unique()
    scrambles = clubs_df['scramble'].unique()

    for block_size in block_sizes:
        for scramble in scrambles:
            filtered_clubs_df = clubs_df[(clubs_df['block_size'] == block_size) & (clubs_df['scramble'] == scramble)]
            filtered_metis_df = metis_df[(metis_df['block_size'] == block_size) & (metis_df['scramble'] == scramble)]
            filtered_original_df = original_df[(original_df['block_size'] == block_size) & (original_df['scramble'] == scramble)]


            if filtered_clubs_df.empty and filtered_metis_df.empty:
                continue

            clubs_best = get_best_results(filtered_clubs_df, ['matrix_name'])
            metis_best = get_best_results(filtered_metis_df, ['matrix_name'])
            original_best = get_best_results(filtered_original_df, ["matrix_name"])

            geometric_mean_ratios(original_best, clubs_best, 'clubs')
            geometric_mean_ratios(original_best, metis_best, 'metis')

            combined_best = pd.concat([original_best,clubs_best, metis_best])

            pivot_df = combined_best.pivot_table(index='matrix_name', columns='algo', values='VBR_nzblocks_count', fill_value=0)
            pivot_df.plot(kind='bar', figsize=(14, 8))

            plt.title(f'Block Size: {block_size}, Scramble: {scramble}')
            plt.xlabel('Matrix Name')
            plt.ylabel('VBR nzblocks count')
            plt.legend(title='Algorithm')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()

            plot_filename = os.path.join(output_dir, f'barplot_bs{block_size}_scramble{scramble}.png')
            plt.savefig(plot_filename)
            plt.close()

# Function to calculate the geometric mean of ratios
def geometric_mean_ratios(original_df, processed_df, algo_name):
    ratios = []
    for matrix_name in processed_df['matrix_name'].unique():
        original_nzblocks = original_df[original_df['matrix_name'] == matrix_name]['VBR_nzblocks_count'].values
        processed_nzblocks = processed_df[processed_df['matrix_name'] == matrix_name]['VBR_nzblocks_count'].values
        if original_nzblocks.size > 0 and processed_nzblocks.size > 0:
            ratio = processed_nzblocks[0] / original_nzblocks[0]
            ratios.append(ratio)
    if ratios:
        geomean = np.exp(np.mean(np.log(ratios)))
        print(f"Geometric mean of {algo_name} ratios: {geomean}")
    else:
        print(f"No valid ratios to calculate for {algo_name}")


# Create bar plots
create_bar_plots(clubs_df, metis_df, original_df, output_plot_dir)







