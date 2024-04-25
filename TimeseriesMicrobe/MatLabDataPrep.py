#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 07:14:29 2022

@author: bkamos & mjgordo
"""
#%%
import numpy as np
import pandas as pd
from scipy.stats import truncnorm, norm, rankdata
import seaborn as sns
from skbio.stats.composition import clr
import matplotlib.pyplot as plt


#%% create a gaussian represenation of the data through rank ordering
def gaussianize(x_in, shuffle = False):
    # rank gaussianize data x
    # rows = samples, columns = features
    x = np.array(x_in)
    x_gauss = np.zeros_like(x)
    if shuffle == True:
        x = np.random.default_rng().permuted(x)
    for feature in range(x.shape[1]):
        x_rank = rankdata(x[:,feature],method='ordinal',nan_policy='omit')
        x_rank = (x_rank-0.5)/np.sum(~np.isnan(x_rank)) #len but nan-robust
        x_gauss[:,feature] = norm.ppf(x_rank)
    return x_gauss


#%% import the data
data_raw = pd.read_csv('/FilteredData.csv')
taxa = pd.read_csv('/taxonomy.csv', sep = ',')
meta = pd.read_csv('/ES_Meta.csv')


# %%
data_raw.columns = data_raw.columns.str.replace('B','')
data = data_raw.set_index('Sample')
data = data.T
data['Day'] = data.index.str.split('_').str[1]
data['Day'] = pd.to_numeric(data['Day'])
data['Rep'] = data.index.str.split('_').str[2]
data['Rep'] = pd.to_numeric(data['Rep'])
data = data.sort_values(by=['Rep', 'Day'])
days = sorted(list(set(data['Day'])))
reps = sorted(list(set(data['Rep'])))


#%% Drop day and rep columns and clr transform the whole frame at once
data_cols = data[['Day','Rep']]
data = data.drop(columns=['Day','Rep'])
data_header = data.columns
samples = data.index
data += 1e-6
data_clr = clr(data)
clr_df = pd.DataFrame(data=data_clr, columns=data_header, index=samples )
clr_df = pd.concat([data_cols, clr_df], axis=1)
clr_df = clr_df.drop(columns=['Day'])
clr_df = clr_df.rename(columns={'Rep': 'Subject'})


#%% Write out the specific frames
data_123 = clr_df[clr_df['Subject'].isin([1,2,3])]
data_234 = clr_df[clr_df['Subject'].isin([2,3,4])]
data_134 = clr_df[clr_df['Subject'].isin([1,3,4])]
data_124 = clr_df[clr_df['Subject'].isin([1,2,4])]
data_one = clr_df[clr_df['Subject'].isin([1])]
data_two = clr_df[clr_df['Subject'].isin([2])]
data_three = clr_df[clr_df['Subject'].isin([3])]
data_four = clr_df[clr_df['Subject'].isin([4])]


#%% isolate the dataframes by rep and gaussinize by rep
grouped = data.groupby('Rep')
rep_dfs = []
for group_name, group_df in grouped:
    # 'group_name' will contain the unique value from the column
    # 'group_df' will contain the DataFrame corresponding to that value
    rep_dfs.append(group_df)

gaussianized_Data = []
day_rep_cols = []
for i in rep_dfs:
    samples = i.index
    data_cols = i[['Day','Rep']]
    i = i.drop(columns=['Day','Rep'])
    isolates = i.columns
    data_ranked = gaussianize(i)
    data = pd.DataFrame(data=data_ranked, columns=isolates, index=samples)
    gaussianized_Data.append(data)
    day_rep_cols.append(data_cols)

#%% combine data that has been gaussinized by batch
concatenated_dfs = [pd.concat([day_rep_cols, gaussianized_Data], axis=1) for day_rep_cols, gaussianized_Data in zip(day_rep_cols, gaussianized_Data)]
concatenated_dfs = pd.concat(concatenated_dfs, axis=0)

#%% Write out gaussanized by batch data
concat_df = concatenated_dfs.drop(columns=['Day'])
concat_df = concat_df.rename(columns={'Rep': 'Subject'})

# %% Write out the specific data frames
data_123 = concat_df[concat_df['Subject'].isin([1,2,3])]
data_234 = concat_df[concat_df['Subject'].isin([2,3,4])]
data_134 = concat_df[concat_df['Subject'].isin([1,3,4])]
data_124 = concat_df[concat_df['Subject'].isin([1,2,4])]
data_one = concat_df[concat_df['Subject'].isin([1])]
data_two = concat_df[concat_df['Subject'].isin([2])]
data_three = concat_df[concat_df['Subject'].isin([3])]
data_four = concat_df[concat_df['Subject'].isin([4])]

#%% get data in a rank order and extract the information of importance for all the data at once, irregardless of batch
samples = data.index
data_cols = data[['Day','Rep']]
data = data.drop(columns=['Day','Rep'])
isolates = data.columns
data_ranked = gaussianize(data)
data = pd.DataFrame(data=data_ranked, columns=isolates, index=samples)

#%% Isolate out only the top 20 most abundant isolates
# Calculate the sum of each column
sums = data.sum()

# Find the column(s) with the highest sum
max_sum_columns = sums.nlargest(20).index

# Create a new DataFrame with only the columns with the highest sum
new_df_max_iso = data[max_sum_columns]


# %% Calculatae the relative abundance of each microbe by row (sample)
row_sums = new_df_max_iso.sum(axis=1)
df_percentages = new_df_max_iso.div(row_sums, axis=0)


# %%
row_sums_check = df_percentages.sum(axis=1)

# %%
data_1 = pd.concat([data_cols, data], axis=1)
data_1 = data_1.drop(columns=['Day'])
data_1 = data_1.rename(columns={'Rep': 'Subject'})

# %%
data_123 = data_1[data_1['Subject'].isin([1,2,3])]
data_234 = data_1[data_1['Subject'].isin([2,3,4])]
data_134 = data_1[data_1['Subject'].isin([1,3,4])]
data_124 = data_1[data_1['Subject'].isin([1,2,4])]
data_one = data_1[data_1['Subject'].isin([1])]
data_two = data_1[data_1['Subject'].isin([2])]
data_three = data_1[data_1['Subject'].isin([3])]
data_four = data_1[data_1['Subject'].isin([4])]


# %%
