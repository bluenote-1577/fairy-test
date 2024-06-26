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

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

dirt = sys.argv[1:]


cm = 1/2.54  # centimeters in inches\n",
    ##Change this to get a bigger figure. \n",
plt.rcParams.update({'font.size': 7})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

# Function to read and label the data
def read_and_label(file):
    # Read the TSV file
    data = pd.read_csv(file, sep='\t')
    
    # Rename the columns to make them concordant if needed
    # This is a generic renaming, adjust according to your actual data
    # For example, if columns in type B files contain additional path information like 'scaffolds.fasta/', you would remove it
    data.columns = [col.replace('QC/reads/', '').replace('.fastq.gz', '').replace('.fastq.gz-var', '_var').replace('scaffolds.fasta/', '') for col in data.columns]
    
    if 'wallen' in file:
        data['Dataset'] = 'Human Gut (n = 10)'
    elif 'sediment' in file:
        data['Dataset'] = 'Sediment (n = 12)'
    elif 'soil' in file or 'olm' in file:
        data['Dataset'] = 'Soil (n = 10)'
    elif 'spmp' in file:
        data['Dataset'] = 'Human Gut - Nanopore Old (n = 7)'
    elif 'sludge' in file:
        data['Dataset'] = 'Sludge - Nanopore R10.4 (n = 5)'
    elif 'pacbio' in file:
        data['Dataset'] = 'Water - PacBio HiFi (n = 4)'
    elif 'biofilm' in file:
        data['Dataset'] = 'Biofilm (n = 8)'


    if 'bwa' in file:
        data['Coverage'] = 'BWA'
    elif 'minimap' in file:
        data['Coverage'] = 'minimap2'
    else:
        data['Coverage'] = 'fairy'

    if 'metabat' in file:
        data['Binner'] = 'MetaBAT2'
    else:
        data['Binner'] = 'VAMB'

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

# Define the directory paths to the type A and type B files

# List all the files in the directories
type_a_files = []
type_b_files = []
for d in dirt:
    for currentpath, folders, files in os.walk(d):
        for file in files:
            f = os.path.join(currentpath, file)
            if 'quality_report' in f:
                if 'fairy' in f:
                    type_b_files.append(f)
                else:
                    type_a_files.append(f)

type_a_files.sort()
type_b_files.sort()
# Read and label each file, then concatenate them into a single DataFrame
print(type_a_files)
type_a_data = [read_and_label(f) for f in type_a_files]
type_b_data = [read_and_label(f) for f in type_b_files]

# Concatenate all data into a single DataFrame
df = pd.concat(type_a_data + type_b_data, ignore_index=True)
print(df.head())
df['High-quality'] = (df['Completeness'] > 90) & (df['Contamination'] < 5)
df['MH-quality'] = (df['Completeness'] > 70) & (df['Contamination'] < 5)
df['Medium-quality'] = (df['Completeness'] > 50) & (df['Contamination'] < 5)

sns.color_palette("light:b", as_cmap=True)
# Create a new DataFrame from these counts
quality_counts = df.groupby(['Sampling', 'Binner', 'Coverage', 'Technology', 'Dataset'])[['High-quality', 'MH-quality', 'Medium-quality']].sum().reset_index()
fig, axs = plt.subplots(4,2, figsize=(12*cm, 12*cm), squeeze=False)

plot_groups = quality_counts.groupby(['Dataset'], sort=False)
for k_p,(title, quality_counts) in enumerate(plot_groups):
    k = k_p % 2
    quality_counts_single = quality_counts[quality_counts['Binner'] == 'VAMB']
    quality_counts_multi = quality_counts[quality_counts['Binner'] == 'MetaBAT2']
    print(quality_counts_single)
    contains_minimap = (quality_counts['Technology'] == 'Nanopore').any()
    if contains_minimap:
        hue_order = ['fairy', 'minimap2']
    else:
        hue_order = ['fairy', 'BWA']

    kf = 2 * (k_p//2)

    width = 0.8
    ql = ['Medium-quality','MH-quality', 'High-quality']
    pals = [sns.color_palette('pastel'), 
            sns.color_palette('muted'), 
            sns.color_palette('dark')]

    if not contains_minimap:
        ss = [0,1,2]
    else:
        ss = [0,2,1]
    pals = [[x[s] for s in ss] for x in pals]

    for i in range(kf, 2 +kf):
        if i%2 == 0:
            df_p = quality_counts_multi
        else:
            df_p = quality_counts_single
        for j in range(len(ql)):
            print(df_p, 'test')
            leg = False
            if j == 1 and i == 0 and k == 1:
                leg = True
            sns.barplot(data = df_p,
                    hue_order = hue_order,
                    x = ql[j],
                    y='Sampling',
                    hue='Coverage',
                    palette = pals[j],
                    legend=False,
                    saturation=1.0,
                    width = width,
                    ax = axs[i][k])
            if i == 0:
                axs[i][k].set_yticklabels(['MetaBAT2\nMulti', 'MetaBAT2\nSingle'])
            else:
                axs[i][k].set_yticklabels(['VAMB\nMulti', 'VAMB\nSingle'])

    axs[0 + kf][k].set_title(title[0])
    axs[1 + kf][k].sharex(axs[kf][k])
    axs[1 + kf][k].set_xlabel("> 90%, 70%, 50% complete\n and <5 % contamination")
    axs[0 + kf][k].set_xlabel("")
#    if k > 0:
#        axs[0][k].set_yticklabels("")
#        axs[1][k].set_yticklabels("")


#lines_labels = [fig.axes]
#lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
#fig.legend(lines, labels,ncol = len(dirt))
muted = sns.color_palette('muted')
import matplotlib.patches as mpatches
labels = ['fairy', 'BWA', 'minimap2']
patches = [mpatches.Patch(color=muted[i], label = labels[i]) for i in range(len(labels))]
#axs[0][1].legend(handles = patches, frameon=False, loc= 'upper center', ncol = len(dirt), bbox_transform = plt.gcf().transFigure, bbox_to_anchor = (0.5, 1.00))
plt.rcParams.update({'font.size': 6.5})
for a in axs:
    for ax in a:
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_ylabel("")
#        for i,cont in enumerate(ax.containers):
#            if i == 0 or i == 1:
#                ax.bar_label(cont, fmt = '%g', label_type = 'edge', padding=1.5)
#            else:
#                ax.bar_label(cont, fmt = '%g', label_type = 'edge', padding=0)
#            
plt.tight_layout()
plt.show()
