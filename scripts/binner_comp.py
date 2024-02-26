import sys
import os
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
        data['Dataset'] = 'Human Gut'
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
    else :
        data['Binner'] = 'MetaBinner'

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
all_dirs = [metabinner_dirs, maxbin_dirs, vamb_dirs, metabat_dirs]
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

quality_counts = df.groupby(['Sampling', 'Binner', 'Coverage', 'Technology', 'Dataset'])[['High-quality', 'MH-quality', 'Medium-quality']].sum().reset_index()

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

        #high_quality_ratio = (fairy_data['High-quality'].values[0] - bwa_data['High-quality'].values[0]) / (bwa_data['High-quality'].values[0])
        #mh_quality_ratio = (fairy_data['MH-quality'].values[0] - bwa_data['MH-quality'].values[0]) / (bwa_data['MH-quality'].values[0])
        #medium_quality_ratio = (fairy_data['Medium-quality'].values[0] - bwa_data['Medium-quality'].values[0]) / (bwa_data['Medium-quality'].values[0])

        #high_quality_ratio = (fairy_data['High-quality'].values[0] - bwa_data['High-quality'].values[0]) 
        #mh_quality_ratio = (fairy_data['MH-quality'].values[0] - bwa_data['MH-quality'].values[0]) 
        #medium_quality_ratio = (fairy_data['Medium-quality'].values[0] - bwa_data['Medium-quality'].values[0]) 


        # Append the ratios along with the group name
        ratios_list.append((*name, high_quality_ratio, mh_quality_ratio, medium_quality_ratio))

# Convert the list to a DataFrame
ratios_df = pd.DataFrame(ratios_list, columns=['Sampling', 'Binner', 'Technology', 'Dataset', 'High-Quality Ratio', 'MH-Quality Ratio', 'Medium-Quality Ratio'])

# Display the DataFrame
print(ratios_df)

#seaborn bar plot, hue by dataset, group by binner
# Set the figure size
plt.figure(figsize=(16/2.54, 6.8/2.54))

# Create the bar plot
ratio = 'Medium-Quality Ratio'
#ratio = 'High-Quality Ratio'
bar = sns.barplot(data=ratios_df, x='Binner', y=ratio, hue='Dataset', ci=None, palette=sns.color_palette('hls'))
#bar.set(ylim=(0.5, 1.2))
plt.legend(frameon=False, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.1))

#print the average of medium-quality ratio for each binner
print(ratios_df.groupby('Binner')[ratio].mean())

# Set the title and labels
#plt.title('High-Quality Ratio')
plt.ylabel("Percentage of MAGs\nrecovered compared to BWA")
#plt.title("Fairy vs BWA: > 50% completeness and < 5% contamination (multi-coverage)")
plt.axhline(100, color='black', linestyle='--', linewidth=1)

#add a grid
plt.grid(axis='y')
#plt.ylim(0, 130)

#logscale
#plt.yscale('log')

#remove upper and right borders
sns.despine()

plt.savefig('figures/fig2-bincomp.svg', format='svg', dpi=300, bbox_inches='tight')
plt.xlabel("")
# Show the plot
plt.show()
