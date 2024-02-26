import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns

folders = ['benchmarks1', 'benchmarks3', 'benchmarks5', 'benchmarks7', 'benchmarks9']


cm = 1/2.54
fs = 7
plt.rcParams.update({'font.size': fs})
#plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})

times = [1,3,5,7,9]
df_all = pd.DataFrame()
for i,folder in enumerate(folders):
    #walk through folder and read every file into a dataframe
    for file in os.listdir(folder):
        if 'bwa_all' in file:
            print(file, 'bwa-all')
            df = pd.read_csv(os.path.join(folder, file),sep='\t')
            df['Coverage tool'] = 'BWA all-to-all'
        elif 'fairy' in file:
            print(file,' fairy')
            df = pd.read_csv(os.path.join(folder, file),sep='\t')
            df['Coverage tool'] = 'fairy all-to-all'
        elif 'bwa_single' in file:
            print(file, 'bwa-single')
            df = pd.read_csv(os.path.join(folder, file),sep='\t')
            df['Coverage tool'] = 'BWA single-cov'

        df['Number of samples'] = times[i]
        df_all = pd.concat([df_all, df], ignore_index=True)


fig,ax = plt.subplots(1,1,figsize=(5.5 * cm, 5.5 * cm))
# stratify by coverage tool and plot the sum of the times as a lineplot with x axis number of samples
df_total_time = df_all.groupby(['Coverage tool', 'Number of samples']).sum().reset_index()
print(df_total_time)
#put markers on the points
#log scale
#ax.set_yscale('log')
sns.lineplot(data=df_total_time, x='Number of samples', y='s', hue='Coverage tool', ax=ax, marker='o', 
        palette=sns.color_palette('muted'), hue_order = ['fairy all-to-all', 'BWA all-to-all', 'BWA single-cov'])
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Wall time (s)\n40 threads')
plt.legend(frameon=False)

# annotate the points, offset it slightly, and cut to 0 decimal places
for i, row in df_total_time.iterrows():
    #don't annotate if it's the first sample time
    if row['Number of samples'] != 9:
        continue
    else:
        ax.annotate(f"{row['s']:.0f}", (row['Number of samples'], row['s']), xytext=(0, 5), textcoords='offset points', ha='center', va='baseline')


#ax.set_title("Multi-sample coverage Pearson R\n between fairy and BWA for contigs")
plt.tight_layout()
plt.savefig('figures/timing.svg')
plt.show()


#plt.show()

    #print(df1.head())
    #print(df2.head())


