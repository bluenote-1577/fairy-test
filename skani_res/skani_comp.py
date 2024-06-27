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
        data['Dataset'] = 'Biofilm'
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

all_dirs = [metabat_dirs]
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

print(df.head())


# skani outputs look like 
#Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
#biofilm/bins/SRR6869036_fairy_all_metabat2_short/SRR6869036_fairy_all_metabat2_short_bin.16.fa	biofilm/bins/SRR6869036_bwa_all_metabat2_short/SRR6869036_fairy_all_metabat2_short_bin.36.fa	100.00	29.16	94.73	NODE_1002_length_31263_cov_13.180947	NODE_1002_length_31263_cov_13.180947
#biofilm/bins/SRR6869035_fairy_all_metabat2_short/SRR6869035_fairy_all_metabat2_short_bin.32.fa	biofilm/bins/SRR6869035_bwa_all_metabat2_short/SRR6869035_fairy_all_metabat2_short_bin.118.fa	100.00	39.60	3.44	NODE_100_length_79714_cov_9.756738	NODE_100_length_79714_cov_9.756738
#biofilm/bins/SRR6869035_fairy_all_metabat2_short/SRR6869035_fairy_all_metabat2_short_bin.88.fa	biofilm/bins/SRR6869035_bwa_all_metabat2_short/SRR6869035_fairy_all_metabat2_short_bin.118.fa	100.00	95.96	87.84	NODE_657_length_40951_cov_8.982590	NODE_100_length_79714_cov_9.756738

# read the skani output
envs = ['wallen_subsamp', 'olm_soil', 'sediment', 'biofilm', 'chicken']
d = {"wallen_subsamp": "Human gut", "olm_soil": "Soil", "sediment": "Sediment", "biofilm": "Biofilm", "chicken": "Chicken caecum"}
skani_df = pd.DataFrame()
for dir in envs:
    skani = pd.read_csv('skani_res/skani_' + dir + '_metabat2.tsv', sep='\t')
    skani['Environment'] = d[dir]
    #take last part of file path and remove extension
    #skani['Ref_name'] = skani['Ref_file'].apply(lambda x: x.split('/')[-1].split('.')[0:-1].join('.'))
    skani['Ref_name'] = skani['Ref_file'].apply(lambda x: '.'.join(x.split('/')[-1].split('.')[0:-1]))
    skani['Query_name'] = skani['Query_file'].apply(lambda x: '.'.join(x.split('/')[-1].split('.')[0:-1]))

    skani['Same_sample'] = skani['Query_name'].apply(lambda x: x.split('_')[0]) == skani['Ref_name'].apply(lambda x: x.split('_')[0])
    #query_sample = skani['Query_name'].apply(lambda x: x.split('_')[0])
    #ref_sample = skani['Ref_name'].apply(lambda x: x.split('_')[0])
    #print(query_sample, ref_sample)
    #if query_sample.equals(ref_sample):
    #    skani['Same_sample'] = True
    #else:
    #    skani['Same_sample'] = False
    skani_df = pd.concat([skani_df, skani])


# keep rows in skani_df with Query_name being medium-quality in df
skani_df = skani_df[skani_df['Same_sample']]
skani_df = skani_df[skani_df['Query_name'].isin(df[df['Medium-quality']]['Name'])]

# for each query_name, keep the row with highest AF and > 99.5% ANI
skani_df = skani_df[skani_df['ANI'] > 99.0]
skani_df = skani_df.loc[skani_df.groupby('Query_name')['Align_fraction_ref'].idxmax()]
# save 
skani_df.to_csv('skani_res/skani_df.tsv', sep='\t', index=False)

# violin plot of AF values for each environment 
fig, ax = plt.subplots(figsize=(14 * 1/2.54,6 * 1/2.54))
ord = ['Human gut', 'Biofilm', 'Sediment', 'Soil', 'Chicken caecum']
#sns.boxplot(x='Environment', y='Align_fraction_ref', data=skani_df, ax=ax)
sns.violinplot(x='Environment', y='Align_fraction_ref', data=skani_df, ax=ax, cut=0, width=0.8, inner=None, linewidth=1.5, hue = 'Environment', palette = 'hls', order = ord, hue_order = ord)

sns.despine()
# label medians above violins, sort environment by order
#medians = skani_df.groupby(['Environment'])['Align_fraction_ref'].median()
medians = []

for env in ord:
    medians.append(skani_df[skani_df['Environment'] == env]['Align_fraction_ref'].median())

for i, median in enumerate(medians):
    ax.text(i, 102, f"{median:.2f}", horizontalalignment='center', color='black')

plt.ylabel('Alignment fraction\n(median % on top)')
plt.title('Fairy vs BWA bin sequence comparison (> 50% complete + MetaBAT2)', fontsize=fs)
plt.xlabel("")

plt.tight_layout()
plt.savefig('figures/fairy_bwa_af_violin.pdf')
plt.show()



