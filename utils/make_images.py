import os
import pandas as pd
import matplotlib.pyplot as plt

# Define the directory containing the files
directory = 'results/results_2024'
output_plot_dir = 'results/results_2024/plots'
os.makedirs(output_plot_dir, exist_ok=True)


# Initialize empty dataframes
clubs_df = pd.DataFrame()
metis_df = pd.DataFrame()

# Column names for the dataframes
clubs_columns = [
    "matrix_name", "rows", "cols", "nnz", "block_size", "scramble", "algo", "mask", "tau", "centroid",
    "VBR_nzcount", "VBR_nzblocks_count", "VBR_average_height", "VBR_longest_row"
]
metis_columns = [
    "matrix_name", "rows", "cols", "nnz", "block_size", "scramble", "algo", "metis_obj", "metis_part",
    "VBR_nzcount", "VBR_nzblocks_count", "VBR_average_height", "VBR_longest_row"
]

# Iterate over files in the directory
for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    
    # Check if the file is a 'clubs' or 'metis' file and read it into the appropriate dataframe
    if filename.startswith("clubs"):
        temp_df = pd.read_csv(filepath, delim_whitespace=True, names=clubs_columns,header=0)
        clubs_df = pd.concat([clubs_df, temp_df], ignore_index=True)
    elif filename.startswith("metis"):
        temp_df = pd.read_csv(filepath, delim_whitespace=True, names=metis_columns,header=0)
        metis_df = pd.concat([metis_df, temp_df], ignore_index=True)

combined_df = pd.concat([clubs_df, metis_df], ignore_index=True)

# Function to get the best result for each matrix
def get_best_results(df, group_cols, value_col='VBR_nzblocks_count'):
    return df.loc[df.groupby(group_cols)[value_col].idxmin()]

# Function to create and save bar plots
def create_bar_plots(clubs_df, metis_df, output_dir):
    block_sizes = clubs_df['block_size'].unique()
    scrambles = clubs_df['scramble'].unique()

    for block_size in block_sizes:
        for scramble in scrambles:
            filtered_clubs_df = clubs_df[(clubs_df['block_size'] == block_size) & (clubs_df['scramble'] == scramble)]
            filtered_metis_df = metis_df[(metis_df['block_size'] == block_size) & (metis_df['scramble'] == scramble)]

            if filtered_clubs_df.empty and filtered_metis_df.empty:
                continue

            clubs_best = get_best_results(filtered_clubs_df, ['matrix_name'])
            metis_best = get_best_results(filtered_metis_df, ['matrix_name'])

            combined_best = pd.concat([clubs_best, metis_best])
            combined_best['algo'] = combined_best['algo'].map(lambda x: 'clubs' if x == 'clubs' else 'metis')

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

# Create bar plots
create_bar_plots(clubs_df, metis_df, output_plot_dir)