import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def process(s):
    return s.split('/')[-1].split(".")[0]


def read_reorder_label(df, var = False):
    df = df[df.columns.drop(list(df.filter(regex='-var')))]
    df = df.rename(columns=process)
    df = df.reindex(sorted(df.columns), axis=1)
    return df


file1s = ['pearson_covs/SRR6869032_bwa_all_short.tsv',
'pearson_covs/SRR19064410_bwa_all_short.tsv',
'pearson_covs/DRR310871_bwa_all_short.tsv',
'pearson_covs/SRR2546421_bwa_all_short.tsv',
'pearson_covs/ERR7625420_minimap_all_long.tsv',
'pearson_covs/ERR7014844_minimap_all_long.tsv']

file2s = ['pearson_covs/SRR6869032_fairy_all_short.tsv',
'pearson_covs/SRR19064410_fairy_all_short.tsv',
'pearson_covs/DRR310871_fairy_all_short.tsv',
'pearson_covs/SRR2546421_fairy_all_short.tsv',
'pearson_covs/ERR7625420_fairy_all_long.tsv',
'pearson_covs/ERR7014844_fairy_all_long.tsv']

datasets = ['Biofilm', 'Human Gut', 'Sediment', 'Soil', 'Sludge (Nanopore)', 'Human Gut (Nanopore)']
plt_df = pd.DataFrame([])

medians = []
for f1,f2,dataset in zip(file1s,file2s, datasets):
    df1 = pd.read_csv(f1, sep='\t')
    df2 = pd.read_csv(f2, sep='\t')
    df1 = read_reorder_label(df1)
    df2 = read_reorder_label(df2)
    print(df1.head())
    # Extract columns with "SRR" in their names
    df1_srr_cols = df1.filter(like='RR')
    df2_srr_cols = df2.filter(like='RR')
    # Transpose the DataFrames to switch rows and columns
    df1_transposed = df1.filter(like='RR').T
    df2_transposed = df2.filter(like='RR').T

    # Calculate Pearson correlation between the transposed DataFrames
    correlation_matrix = df1_transposed.corrwith(df2_transposed)
    # Create a new DataFrame to store the correlation matrix
    correlation_df = pd.DataFrame({'Correlation': correlation_matrix})

    # Add a new column 'Dataset' with 'DATASET1' values
    correlation_df['Dataset'] = dataset

    # Reset the index to match the row index of the original DataFrames
    correlation_df.reset_index(drop=True, inplace=True)

    # Append the correlation DataFrame to the existing DataFrame
    plt_df = pd.concat([plt_df, correlation_df], axis=0)

    #print the median correlation
    medians.append(correlation_df['Correlation'].median())

cm = 1/2.54
fs = 7
plt.rcParams.update({'font.size': fs})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

fig,ax = plt.subplots(1,1,figsize=(7.5 * cm, 5.5 * cm))
plt_df = plt_df.dropna()
sns.violinplot(x = 'Dataset', y='Correlation', data = plt_df, ax = ax, hue = 'Dataset', inner=None, palette = sns.color_palette('hls'))
ax.set_ylabel('Pearson R between contig coverage')
ax.spines[['right', 'top']].set_visible(False)
ax.set_title("BWA vs fairy, multi-sample coverage comparison", fontsize=fs)
ax.set_xlabel('')
#annotate the medians above the violins
for i, median in enumerate(medians):
    ax.text(i, 0.08 + median, f"{median:.3f}", horizontalalignment='center', color='black')
plt.tight_layout()
plt.savefig('figures/pearson.svg')



plt.show()

    #print(df1.head())
    #print(df2.head())


