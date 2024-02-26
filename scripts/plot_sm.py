import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
cmap = sns.color_palette('muted')

# Aggregated DataFrames for each type
aggregate_bwa_mags = pd.DataFrame()
aggregate_single_mags = pd.DataFrame()
aggregate_mags = pd.DataFrame()
aggregate_bwa_single_mags = pd.DataFrame()

# Function to read and process the data
def process_data(folder):
    file_path = os.path.join(folder, 'quality_report.tsv')
    binner = 'metabat'
    if 'vamb' in folder:
        binner = 'vamb'
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        # Calculate the metric Completeness - 4*Contamination
        df['Metric'] = df['Completeness'] - 5 * df['Contamination']
        df['Binner'] = binner
        return df
    else:
        return None

# Function to plot the ranked metric
def plot_ranked_metric(df_aggregate, folder_type,i):
    df_aggregate = df_aggregate[df_aggregate['Completeness'] > 50]
    df_aggregate = df_aggregate[df_aggregate['Metric'] > 0]
    print(len(df_aggregate[df_aggregate['Metric'] > 50]))
    # Sort the dataframe based on the metric
    df_sorted = df_aggregate.sort_values(by='Metric', ascending=False)
    # Create the rank plot
    l = folder_type

    df_vamb = df_sorted[df_sorted['Binner'] == 'vamb']
    df_mb = df_sorted[df_sorted['Binner'] == 'metabat']
    plt.plot(np.arange(len(df_vamb['Metric'])), df_vamb['Metric'], marker='o', linestyle='--', label = l, ms = 3, c = cmap[i])
    plt.plot(np.arange(len(df_mb['Metric'])), df_mb['Metric'], marker='o', linestyle='-', label = l, ms = 3, c = cmap[i])
    plt.title(f'Aggregated Ranked Metric Plot for {folder_type}')
    plt.xlabel('Rank')
    plt.ylabel('Completeness - 5 * Contamination')
    plt.grid(True)

# Iterate through each folder in the base directory
bwa_count = 0
fairy_count = 0
for folder_name in sys.argv[1:]:
    if 'bwa_single' in folder_name or 'minimap_single' in folder_name:  # This will match non-BWA mags folders
        df = process_data(folder_name)
        if df is not None:
            aggregate_bwa_single_mags = pd.concat([aggregate_bwa_single_mags, df])
    elif 'bwa_all' in folder_name or 'minimap_all' in folder_name:
        df = process_data(folder_name)
        if df is not None:
            bwa_count += 1
            aggregate_bwa_mags = pd.concat([aggregate_bwa_mags, df])
    elif 'fairy_single' in folder_name:  # This will match non-BWA mags folders
        df = process_data(folder_name)
        if df is not None:
            aggregate_single_mags = pd.concat([aggregate_single_mags, df])
    elif 'fairy_all' in folder_name:
        df = process_data(folder_name)
        if df is not None:
            fairy_count += 1
            aggregate_mags = pd.concat([aggregate_mags, df])

    

# Plot the aggregated metrics
plt.figure(figsize=(10,6))
if not aggregate_bwa_mags.empty:
    plot_ranked_metric(aggregate_bwa_mags, 'bwa-all', 0)

if not aggregate_mags.empty:
    plot_ranked_metric(aggregate_mags, 'fairy-all', 1)

if not aggregate_bwa_single_mags.empty:
    plot_ranked_metric(aggregate_bwa_single_mags, 'bwa-single', 2)

if not aggregate_single_mags.empty:
    plot_ranked_metric(aggregate_single_mags, 'fairy-single', 3)

print('bwa-all', bwa_count)
print('fairy-all', fairy_count)
#plt.ylim([0,100])
#plt.ylim([0,100])
plt.legend()
plt.axhline(50, ls = '--')
plt.show()
