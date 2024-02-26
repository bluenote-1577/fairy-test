import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Read the 'batchfile.tsv' and create a dictionary for ID to filename mapping
def get_df(batchfile, resfile):
    batchfile_df = pd.read_csv(batchfile, sep='\t', header=None, names=['filename', 'ID'])
    batch_dict = dict(zip(batchfile_df['ID'], batchfile_df['filename']))

    # Step 2: Read the file 'x' to get the list of IDs
    with open(resfile, 'r') as file:
        ids = file.read().splitlines()

    # Step 3: Find corresponding filenames for each ID in 'x'
    filenames = [batch_dict[id] for id in ids if id in batch_dict]

    # Step 4: Read and extract data from each corresponding 'quality_report.tsv'
    rows = []
    for filename in filenames:
        # Assuming the structure of the filename to construct the path for 'quality_report.tsv'
        quality_report_path = f"../sediment_new/{filename.split('/')[2]}_cm2/quality_report.tsv"
        quality_df = pd.read_csv(quality_report_path, sep='\t')
        #print(filename, quality_df, quality_report_path)
        # Extract the row that matches the filename
        row = quality_df[quality_df['Name'] == filename.split('/')[3][0:-3]].iloc[0]
        single = ""
        if 'single' in batchfile:
            single = ' single'
        if 'bwa' in filename:
            print(filename)
            row['Coverage']  = 'BWA' + single
        else:
            print(filename,'fairy')
            row['Coverage']  = 'fairy' + single
        rows.append(row)

    # Step 5: Create a final DataFrame
    final_df = pd.DataFrame(rows)
    return final_df
df1 = get_df('./batchfile.tsv', './asgard_id_all.tsv')
df2 = get_df('./batchfile_single.tsv', './asgard_id_single.tsv')
df = pd.concat([df1,df2])
df.to_csv('asgard_df.tsv',sep='\t',index=False,header=True)
df['Metric'] = df['Completeness'] - 5 * df['Contamination']
df['Rank'] = df.groupby('Coverage')['Metric'].rank(method='min', ascending=False)
sns.lineplot(x='Rank', y='Metric', hue = 'Coverage', data=df, markers=True, dashes=False, style="Coverage")
plt.ylim([0,100])
#final_df['Contamination'] = final_df['Contamination']/100
#final_df['Completeness'] = final_df['Completeness']/100
#sns.lmplot(data = df,x='Contamination', y = 'Completeness', hue='Coverage', lowess=True)
#sns.lmplot(data = df2, x='Contamination', y = 'Completeness', hue='Coverage', lowess=True)
#sns.scatterplot(data = final_df,x='Contamination', y = 'Completeness', hue='Coverage')
plt.show()

#print the number of rows with completeness > 90 and contamination < 5 for each Coverage and each df
fairy_mult = df1[df1['Coverage'] == 'fairy']
bwa_mult = df1[df1['Coverage'] == 'BWA']
print(fairy_mult)
print(bwa_mult)
print(fairy_mult[(fairy_mult['Completeness'] > 50) & (fairy_mult['Contamination'] < 5)].shape[0])
print(bwa_mult[(bwa_mult['Completeness'] > 50) & (bwa_mult['Contamination'] < 5)].shape[0])


#print(df_1[(df_1['Completeness'] > 90) & (df_1['Contamination'] < 5)].shape[0])
