import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Aggregated DataFrames for each type
aggregate_bwa_mags = pd.DataFrame()
aggregate_single_mags = pd.DataFrame()
aggregate_mags = pd.DataFrame()
aggregate_bc = pd.DataFrame()

# Function to read and process the data
def process_data(folder):
    file_path = os.path.join(folder, 'quality_report.tsv')
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        # Calculate the metric Completeness - 4*Contamination
        df['Metric'] = df['Completeness'] - 5 * df['Contamination']
        return df
    else:
        return None

# Function to plot the ranked metric
def plot_ranked_metric(df_aggregate, folder_type):
    df_aggregate = df_aggregate[df_aggregate['Completeness'] > 50]
    print(len(df_aggregate[df_aggregate['Metric'] > 50]))
    # Sort the dataframe based on the metric
    df_sorted = df_aggregate.sort_values(by='Metric', ascending=False)
    # Create the rank plot
    l = folder_type

    plt.plot(np.arange(len(df_sorted['Metric'])), df_sorted['Metric'], marker='o', linestyle='-', label = l, ms = 3)
    plt.title(f'Aggregated Ranked Metric Plot for {folder_type}')
    plt.xlabel('Rank')
    plt.ylabel('Completeness - 4*Contamination')
    plt.grid(True)

# Iterate through each folder in the base directory
for folder_name in sys.argv[1:]:
    if 'Single' in folder_name:  # This will match non-BWA mags folders
        df = process_data(folder_name)
        if df is not None:
            aggregate_single_mags = pd.concat([aggregate_single_mags, df])
    elif 'bwa' in folder_name or 'minimap' in folder_name:
        df = process_data(folder_name)
        if df is not None:
            aggregate_bwa_mags = pd.concat([aggregate_bwa_mags, df])
    elif 'bc_' in folder_name or 'fairy' in folder_name:  # This will match non-BWA mags folders
        df = process_data(folder_name)
        if df is not None:
            aggregate_bc = pd.concat([aggregate_bc, df])
    elif '_mags' in folder_name:
        df = process_data(folder_name)
        if df is not None:
            aggregate_mags = pd.concat([aggregate_mags, df])

    

# Plot the aggregated metrics
plt.figure(figsize=(10,6))
if not aggregate_bwa_mags.empty:
    plot_ranked_metric(aggregate_bwa_mags, 'bwa')

if not aggregate_mags.empty:
    plot_ranked_metric(aggregate_mags, 'sylph')

if not aggregate_single_mags.empty:
    plot_ranked_metric(aggregate_single_mags, 'single')

if not aggregate_bc.empty:
    plot_ranked_metric(aggregate_bc, 'bc')

plt.ylim([0,100])
plt.legend()
plt.axhline(50, ls = '--')
plt.show()
