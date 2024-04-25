#%%
from pycm import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import truncnorm


#%%
### import the data
data_raw = pd.read_csv('/FilteredData.csv')
taxa = pd.read_csv('/taxonomy.csv', sep = ',')
meta = pd.read_csv('/ES_Meta.csv')


#%% Remove the B value from the sample
data_raw.columns = data_raw.columns.str.replace('B','')
#%% Get data in a format that is suitable for creating random noise
# Num noise columns
num_noise_columns = 40

#%% Extract unique IDs that represent each time point
unique_ids = data_raw.columns.str.split('_').str[:2].str.join('_').unique()
unique_ids = list(unique_ids)
unique_ids = list(map(lambda x: x + '_', unique_ids))
unique_ids.pop(0)


#%% Extract dataframes that represent each timepoint (each DF is the batch rep of the time)
dfs1 = []
for unique_id in unique_ids:
    dfs = pd.DataFrame()
    cols_with_id = [col for col in data_raw.columns if col.startswith(unique_id)]
    dfs[cols_with_id] = data_raw[cols_with_id]
    dfs1.append(dfs)



#%% Creating random gaussian noise for each a left out batch and the remaining batches - seperatly
for i in range(len(dfs1)):
    std = dfs1[i].std(axis=1)
    std[std == 0] = 0.001
    for j in dfs1[i].columns:
        if j.endswith('1') or j.endswith('2') or j.endswith('3') or j.endswith('4'):
            print(dfs1[i][j].name)
            mean = dfs1[i][j]
            for k in range(10):
                noise_col_name = f'{dfs1[i][j].name}_{k+5}'
                dfs1[i][noise_col_name] = truncnorm.rvs((-mean)/std, np.inf, loc=mean, scale= std, size=len(dfs1[i]))

            cols_to_keep = dfs1[i].drop(j, axis =1)
            mean_remainder = cols_to_keep.mean(axis=1)
            std_remainder = cols_to_keep.std(axis=1)
            for q in range(30):
                noise_col_name = f'{dfs1[i][j].name}_{q+5}_R'
                dfs1[i][noise_col_name] = truncnorm.rvs((-mean_remainder)/std_remainder, np.inf, loc=mean_remainder, scale=std_remainder, size=len(dfs1[i]))



#%%Combine all DFs and set the axis as the isolates
combined_noise_dfs = pd.concat(dfs1, axis=1)
isolates = data_raw.iloc[:,0]
combined_noise_dfs = combined_noise_dfs.set_index(isolates)
combined_noise_dfs = combined_noise_dfs.T


#%% add a sample column that is a number - remove day column and move sample column to the furthest left column
split_index = combined_noise_dfs.index.str.split('_')
combined_noise_dfs['Day'] = split_index.str[1]
combined_noise_dfs['Day'] = pd.to_numeric(combined_noise_dfs['Day'])
combined_noise_dfs['Rep'] = split_index.str[2:].str.join('_')
combined_noise_dfs = combined_noise_dfs.sort_values(by=['Rep', 'Day'])

#%% Reorder the dataframe columns
columns = combined_noise_dfs.columns.tolist()

# Move the last column to the first position
columns = [columns[-1]] + columns[:-1]

# Reorder the DataFrame with the new column order
combined_noise_dfs = combined_noise_dfs[columns]
reps = combined_noise_dfs['Rep']
combined_noise_dfs = combined_noise_dfs.drop(columns=['Day', 'Rep'])


#%% Calculate the percentage of each microbe as a part of the total
row_sums = combined_noise_dfs.sum(axis=1)
df_percentages = combined_noise_dfs.div(row_sums, axis=0)
row_sums_check = df_percentages.sum(axis=1)
data_1 = pd.concat([reps, df_percentages], axis=1)