import sys
import os
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import re
from collections import defaultdict
import seaborn as sns
from dataclasses import dataclass
from scipy import stats
import pandas as pd

fs = 7
plt.rcParams.update({'font.size': fs})
plt.rcParams.update({'font.family':'arial'})
barplot = False
scatter = False
cross_comp = True
mq = True

# Function to read and label the data
def read_and_label(file):
    # Read the TSV file
    data = pd.read_csv(file, sep='\t')
    #if 'single' in file and 'fairy' in file:
    #    return pd.DataFrame([])
    # Rename the columns to make them concordant if needed
    # This is a generic renaming, adjust according to your actual data
    # For example, if columns in type B files contain additional path information like 'scaffolds.fasta/', you would remove it
    data.columns = [col.replace('QC/reads/', '').replace('.fastq.gz', '').replace('.fastq.gz-var', '_var').replace('scaffolds.fasta/', '') for col in data.columns]
    
    if 'wallen' in file:
        data['Dataset'] = 'Human gut'
    elif 'sediment' in file:
        data['Dataset'] = 'Sediment'
    elif 'soil' in file or 'olm' in file:
        data['Dataset'] = 'Soil'
    elif 'spmp' in file:
        data['Dataset'] = 'Human Gut - Nanopore\n(n = 7)'
    elif 'sludge' in file:
        data['Dataset'] = 'Sludge - Nanopore\n(n = 5)'
    elif 'pacbio' in file:
        data['Dataset'] = 'Water - PacBio HiFi\n(n = 4)'
    elif 'biofilm' in file:
        data['Dataset'] = 'Marine biofilm'
    elif 'chicken' in file:
        data['Dataset'] = 'Chicken caecum'

    if 'bwa' in file:
        data['Coverage'] = 'BWA'
    elif 'minimap' in file:
        data['Coverage'] = 'minimap2'
    elif 'fairy' in file:
        data['Coverage'] = 'fairy'

    if 'metabat' in file:
        data['Binner'] = 'MetaBAT2'
    elif 'vamb' in file:
        data['Binner'] = 'VAMB'
    elif 'maxbin' in file:
        data['Binner'] = 'MaxBin2*' 
    elif 'metabinner' in file:
        data['Binner'] = 'MetaBinner'
    elif 'semibin' in file:
        data['Binner'] = 'SemiBin2'

    if 'single' in file:
        data['Sampling'] = 'Single-coverage'
    else:
        data['Sampling'] = 'Multi-coverage'

    if 'short' in file:
        data['Technology'] = 'Illumina'
    else:
        data['Technology'] = 'Nanopore'

    # Add a new column to label the data
    return data


metabinner_dirs = ['wallen_metabinner','soil_metabinner','sediment_metabinner','biofilm_metabinner','chicken_metabinner']
maxbin_dirs = ['wallen_maxbin','soil_maxbin','sediment_maxbin','biofilm_maxbin']
vamb_dirs = ['wallen_vamb','soil_vamb','sediment_vamb','biofilm_vamb','chicken_vamb']
metabat_dirs = ['wallen_new','soil_new','sediment_new','biofilm_new','chicken_new']
semibin_dirs = ['revision1/wallen_rev1_checkm2_results', 'revision1/soil_rev1_checkm2_results', 'revision1/sediment_rev1_checkm2_results', 'revision1/biofilm_rev1_checkm2_results', 'revision1/chicken_rev1_checkm2_results']

all_dirs = [metabinner_dirs, maxbin_dirs, vamb_dirs, metabat_dirs, semibin_dirs]
all_df = []

for d in all_dirs:
    # read each file in the directory
    # walk the directory
    for dir in d:
        for root, dirs, files in os.walk(dir):
            for file in files:
                if file.endswith('quality_report.tsv'):
                    # read all files and put into a large dataframe
                    data = read_and_label(os.path.join(root, file))
                    all_df.append(data)
df = pd.concat(all_df)

#stratify by Coverage, and count how many entries have Contamination < 5 and Completeness > 90
df['High-quality'] = (df['Contamination'] < 5) & (df['Completeness'] >= 90)
df['MH-quality'] = (df['Contamination'] < 5) & (df['Completeness'] >= 70) #& (df['Completeness'] < 90)
df['Medium-quality'] = (df['Contamination'] < 5) & (df['Completeness'] >= 50)# & (df['Completeness'] < 70)
df = df[df['Sampling'] == 'Multi-coverage']

# divide by the number of samples

quality_counts = df.groupby(['Sampling', 'Binner', 'Coverage', 'Technology', 'Dataset'])[['High-quality', 'MH-quality', 'Medium-quality']].sum().reset_index()
#quality_counts = df.groupby(['Sampling', 'Binner', 'Coverage', 'Technology', 'Dataset'])[['High-quality', 'MH-quality', 'Medium-quality']].mean().reset_index()

sample_sizes = {"Chicken caecum": 24, "Human gut": 10, "Marine biofilm": 8, "Sediment": 12, "Soil": 10}
#divide by sample sizes
quality_counts['High-quality'] = quality_counts['High-quality'] / quality_counts['Dataset'].map(sample_sizes)
quality_counts['MH-quality'] = quality_counts['MH-quality'] / quality_counts['Dataset'].map(sample_sizes)
quality_counts['Medium-quality'] = quality_counts['Medium-quality'] / quality_counts['Dataset'].map(sample_sizes)

df = quality_counts
print(df)

#print df to csv
df.to_csv('quality_counts.csv')


# Group by the relevant columns
grouped = df.groupby(['Sampling', 'Binner', 'Technology', 'Dataset'])

# Initialize a list to store the results
ratios_list = []

# Iterate through each group
for name, group in grouped:
    #exit()
    # Check if both 'BWA' and 'fairy' are in this group
    if 'BWA' in group['Coverage'].values and 'fairy' in group['Coverage'].values:
        # Extract BWA and fairy data
        bwa_data = group[group['Coverage'] == 'BWA']
        fairy_data = group[group['Coverage'] == 'fairy']

        # Calculate ratios for each quality metric
        high_quality_ratio   = 100* (fairy_data['High-quality'].values[0]) / (bwa_data['High-quality'].values[0])
        mh_quality_ratio     = 100* (fairy_data['MH-quality'].values[0]) / (bwa_data['MH-quality'].values[0])
        medium_quality_ratio = 100* (fairy_data['Medium-quality'].values[0]) / (bwa_data['Medium-quality'].values[0])

        fairy_mq = fairy_data['Medium-quality'].values[0]
        bwa_mq = bwa_data['Medium-quality'].values[0]

        fairy_hq = fairy_data['High-quality'].values[0]
        bwa_hq = bwa_data['High-quality'].values[0]


        #high_quality_ratio = (fairy_data['High-quality'].values[0] - bwa_data['High-quality'].values[0]) / (bwa_data['High-quality'].values[0])
        #mh_quality_ratio = (fairy_data['MH-quality'].values[0] - bwa_data['MH-quality'].values[0]) / (bwa_data['MH-quality'].values[0])
        #medium_quality_ratio = (fairy_data['Medium-quality'].values[0] - bwa_data['Medium-quality'].values[0]) / (bwa_data['Medium-quality'].values[0])

        #high_quality_ratio = (fairy_data['High-quality'].values[0] - bwa_data['High-quality'].values[0]) 
        #mh_quality_ratio = (fairy_data['MH-quality'].values[0] - bwa_data['MH-quality'].values[0]) 
        #medium_quality_ratio = (fairy_data['Medium-quality'].values[0] - bwa_data['Medium-quality'].values[0]) 


        # Append the ratios along with the group name
        ratios_list.append((*name, high_quality_ratio, mh_quality_ratio, medium_quality_ratio, fairy_mq, bwa_mq, fairy_hq, bwa_hq))

# Convert the list to a DataFrame
ratios_df = pd.DataFrame(ratios_list, columns=['Sampling', 'Binner', 'Technology', 'Dataset', 'High-Quality Ratio', 'MH-Quality Ratio', 'Medium-Quality Ratio', 'Fairy MQ', 'BWA MQ', 'Fairy HQ', 'BWA HQ'])

dataset_order = ['Human gut', 'Marine biofilm', 'Sediment', 'Soil', 'Chicken caecum']
#### Scatter plots

if scatter:

    #plt.figure(1,1,figsize=(16/2.54, 6.8/2.54))
    # Set the figure size

    #subplot with 4 columns
    fig, axs = plt.subplots(1, 4, figsize=(16/2.54, 5/2.54), sharex=True, sharey=True)

    #plt.figure(figsize=(16/2.54, 6.8/2.54))

    # Create a scatter plot for Fairy MQ vs BWA MQ
    #maker size 12
    #scatter = sns.scatterplot(data=ratios_df, x='BWA MQ', y='Fairy MQ', hue='Dataset', style='Binner', palette=sns.color_palette('hls'), s=100,alpha=0.8)
    #iterate over all binners and plot on axs
    for i, binner in enumerate(ratios_df['Binner'].unique()):
        if binner  == 'MaxBin2*':
            pal = sns.color_palette('hls')[0:4]
        else:
            pal = sns.color_palette('hls')
        if mq:
            scatter = sns.scatterplot(data=ratios_df[ratios_df['Binner'] == binner], x='BWA MQ', y='Fairy MQ', hue='Dataset', palette=pal, s=60,alpha=0.8, ax=axs[i], hue_order=dataset_order)
            axs[i].plot([np.min(ratios_df['BWA MQ']), np.max(ratios_df['BWA MQ'])], [np.min(ratios_df['BWA MQ']), np.max(ratios_df['BWA MQ'])], 'k--', linewidth=1)
        else:
            scatter = sns.scatterplot(data=ratios_df[ratios_df['Binner'] == binner], x='BWA HQ', y='Fairy HQ', hue='Dataset', palette=pal,alpha=0.8, ax=axs[i], hue_order=dataset_order, s=60)
            axs[i].plot([np.min(ratios_df['BWA HQ']), np.max(ratios_df['BWA HQ'])], [np.min(ratios_df['BWA HQ']), np.max(ratios_df['BWA HQ'])], 'k--', linewidth=1)
    #loglog
    for i,ax in enumerate(axs):
        ax.set_xscale('log')
        ax.set_yscale('log')
        if mq:
            ax.set_xlabel("BWA # >50% complete MAGs\n(sample average)")
            ax.set_ylabel("fairy # >50% complete MAGs\n(sample average)")
        else:
            ax.set_xlabel("BWA # >90% complete MAGs\n(sample average)")
            ax.set_ylabel("fairy # >90% complete MAGs\n(sample average)")
        #ax.plot([0.1, 100], [0.1, 100], 'k--', linewidth=1)
        #ax.plot([1, 4e3], [1*1.5, 4e3*1.5], 'k--', linewidth=1)
        ax.legend().remove()
    #    else:
    #        ax.legend(frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.2))
        percent = 50 if mq else 90
        ax.set_title(ratios_df['Binner'].unique()[i] + " (>{percent}% complete)".format(percent=percent)
                     , fontsize=fs)
        sns.despine(ax=ax)
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.plot([0, 4e3], [0, 4e3], 'k--')
    #legend at top
    #plt.gca().add_artist(plt.legend(frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.2)))
    # legend at the top, colours for the dataset, and square for first scatter, circle for second
   # plt.legend(frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.2))
   # Define custom legend handles
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=5, label='Black Circle'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=5, label='Black Square')
    ]
    #grid 
    #plt.grid(axis='both')

    # Create the legend
    #plt.legend(handles=legend_handles,ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.2))

    plt.tight_layout()
    if mq:
        plt.savefig('figures/loglog-mq.svg', format='svg', dpi=300, bbox_inches='tight')
    else:
        plt.savefig('figures/loglog-hq.svg', format='svg', dpi=300, bbox_inches='tight')

    plt.show()

#### Bar plots

if barplot:

    # Display the DataFrame
    print(ratios_df)

    #seaborn bar plot, hue by dataset, group by binner
    # Set the figure size
    plt.figure(figsize=(16/2.54, 6.8/2.54))

    # Create the bar plot
    ratio = 'Medium-Quality Ratio'
    #ratio = 'High-Quality Ratio'
    bar = sns.barplot(data=ratios_df, x='Binner', y=ratio, hue='Dataset', ci=None, palette=sns.color_palette('hls'))
    ratios_df['MQ Difference'] = ratios_df['Fairy MQ'] - ratios_df['BWA MQ']

    # Retrieve the number of unique 'Binner' and 'Dataset' values
    n_binner = 5
    n_dataset = 5

    print(n_binner, n_dataset)

    # Iterate through each bar and annotate
#    print(bar.patches.__len__())
#    for i, p in enumerate(bar.patches):
#        if i >= 19:
#            break
#        if i  >= 16:
#            i+=1
#        # The bar index within the group
#        bar_idx_within_group = i // n_binner
#        # The group index
#        group_idx = i % n_binner
#
#        print(bar_idx_within_group,group_idx)
#
#        # Retrieve the corresponding row from the dataframe
#        row = ratios_df[(ratios_df['Binner'] == ratios_df['Binner'].unique()[group_idx]) & 
#                        (ratios_df['Dataset'] == ratios_df['Dataset'].unique()[bar_idx_within_group])]
#
#
#        if row.empty:
#            continue
#        row = row.iloc[0]
#        print(row)
#        mq_diff = row['MQ Difference']
#
#        # If we have a value to annotate with, do it
#        if mq_diff is not None:
#            # Format the 'MQ Difference' with 2 decimal places
#            if mq_diff < 0:
#                text = f"{mq_diff:.0f}"
#            else:
#                text = f"+{mq_diff:.0f}"
#
#            # Set the position of the text
#            x = p.get_x() + p.get_width() / 2
#            y = p.get_height()
#
#            # Add the text annotation on the top of the bar
#            plt.text(x, y, text, ha='center', va='bottom')
#
    #bar.set(ylim=(0.5, 1.2))

    #print the average of medium-quality ratio for each binner
    print(ratios_df.groupby('Binner')[ratio].mean())

    # Set the title and labels
    #plt.title('High-Quality Ratio')
    plt.ylabel("Percentage of BWA's MAGs recovered\n(> 50% complete, < 5% contaminated)")
    #plt.title("Fairy vs BWA: > 50% completeness and < 5% contamination (multi-coverage)")
    plt.axhline(100, color='black', linestyle='--', linewidth=1)

    #add a grid
    plt.grid(axis='y')
    #plt.ylim(0, 130)

    #logscale
    #plt.yscale('log')

    #remove upper and right borders
    sns.despine()

    plt.legend(frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.1))
    plt.xlabel("")
    plt.savefig('figures/fig2-bincomp.svg', format='svg', dpi=300, bbox_inches='tight')
    # Show the plot
    plt.show()

if cross_comp:
    pals = [sns.color_palette('pastel'), 
            sns.color_palette('muted'), 
            sns.color_palette('dark')]
    quals = ['Medium-quality', 'MH-quality', 'High-quality']

    #create a barplot for each dataset
    #dual bar plot with fairy vs bwa side by side, stratified by binner
    # Set the figure size
    fig, axs = plt.subplots(2, 3, figsize=(16/2.54, 12/2.54))

    # Create the bar plot
    for j in range(3):
        for i, dataset in enumerate(ratios_df['Dataset'].unique()):
            # Create a subplot for each dataset
            dataset_data = df[df['Dataset'] == dataset]

            leg = False
            if i == 1 and j == 1:
                leg = True

            # Create the bar plot for the current dataset
            bar = sns.barplot(data=dataset_data, x='Binner', y=quals[j], hue = 'Coverage', ci=None, palette=pals[j], alpha=0.8, ax = axs[i//3, i%3], legend=leg, hue_order = ['fairy', 'BWA'], saturation=1) 

            #rename legend to "BWA (50, 70 90% complete)"
            handles, labels = axs[i//3, i%3].get_legend_handles_labels()


            
            axs[i//3, i%3].set_ylabel("Number of MAGs\n(sample average)")
            #despine 
            sns.despine(ax=axs[i//3, i%3])

            #remove xlabel
            axs[i//3, i%3].set_xlabel('')

            #rotate xticks
            axs[i//3, i%3].tick_params(axis='x', rotation=45)

            # Set the title and labels
            axs[i//3, i%3].set_title(dataset, fontsize=fs)
            #no legend unless first one
            if i == 1 and j == 1:
                axs[i//3, i%3].legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.35))

    # remove the last subplot
    fig.delaxes(axs[1,2])

            # Set the legend
            
    # sup tittle
    fig.suptitle("Multi-sample short-read binning \n(> 90, 70, and 50% complete and < 5% contaminated)", fontsize=fs)
    plt.tight_layout()
    plt.savefig('figures/binner_lmh.svg', format='svg', dpi=300, bbox_inches='tight')
    plt.show()


