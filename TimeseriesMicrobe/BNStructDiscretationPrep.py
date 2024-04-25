#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 07:14:29 2022

@author: bkamos & mjgordo
"""
#%%
from pycm import *
from skbio.stats.composition import clr
from sklearn.cluster import AgglomerativeClustering, DBSCAN, KMeans
import numpy as np
import pandas as pd
from operator import itemgetter
import csv
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns
from scipy.stats import truncnorm


#%% Define functions of chosing
## The below function has traditionally be run with pd.cut (priarily studied) where it just takes the data range and splits it in half but it would be intersted into check qcut and take the median value and the seperator for bins
def quantile_cut(data):
    num_bins = 2  # Number of bins
    bin_values = []
    for column in data.columns:
        data[f'{column}_bins'] = pd.qcut(data[column], q=num_bins, labels=False)
        # print(data)
        # data[f'{column}_bin_range'] = pd.qcut(data[column], q=num_bins, labels=False, retbins=True)
        bins = pd.qcut(data[column], q=num_bins, labels=False, retbins=True)
        bin_values.append(bins)

    return data, bin_values

def binningZeroOne(data, columnHeadersClean, numberOfBins, listOfLabels, measure):
    
    discretizedDF1 = pd.DataFrame()
    bins1 = []

    for i in range(len(columnHeadersClean)):

        if measure == 'mean':
            discretizedDF1[columnHeadersClean[i]] = pd.cut(x=data[columnHeadersClean[i]], bins=numberOfBins[columnHeadersClean[i]], precision = 3, duplicates = 'drop', labels = listOfLabels, include_lowest=True) # can add numberOfBins[columnHeadersClean[i]] for mean value
                    
            binRange = pd.cut(x=data[columnHeadersClean[i]], bins=numberOfBins[columnHeadersClean[i]], precision = 3, duplicates = 'drop', labels = listOfLabels, retbins=True, include_lowest=True)
        
        else:
            discretizedDF1[columnHeadersClean[i]] = pd.cut(x=data[columnHeadersClean[i]], bins=numberOfBins, precision = 3, duplicates = 'drop', labels = listOfLabels, include_lowest=True) # can add numberOfBins[columnHeadersClean[i]] for mean value
                    
            binRange = pd.cut(x=data[columnHeadersClean[i]], bins=numberOfBins, precision = 3, duplicates = 'drop', labels = listOfLabels, retbins=True, include_lowest=True)
        

        bins1.append(binRange)

    return discretizedDF1, bins1

def sort_tuples(tuples):
    return sorted(tuples, key=itemgetter(1))


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


# #%% Create the random gausian noise for each dataframe based on the mean and std of each isolate - orignal form of the algorithim
# for i in range(len(dfs1)):
#     print(dfs1[i].columns[0])
#     mean = dfs1[i].mean(axis=1)
#     std = dfs1[i].std(axis=1)
#     std[std == 0] = 0.001
#     for j in range(num_noise_columns):
#         noise_col_name = f'{unique_ids[i]}{j+5}'
#         dfs1[i][noise_col_name] = truncnorm.rvs((-mean)/std, np.inf, loc=mean, scale= std, size=len(dfs1[i]))

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


#%% Calculate the threshold based on the new data
flat_dfs_for_thres = combined_noise_dfs.to_numpy().flatten()
threshold = np.percentile(flat_dfs_for_thres, 95)


#%% Convert the orignal noised frame into 0/1 values
df_binary = combined_noise_dfs >= threshold
df_binary = df_binary.astype(int)


#%% Drop the columns that don't change
unique_counts = df_binary.nunique()
cols_to_drop = unique_counts[unique_counts == 1].index
bin_filt = df_binary.drop(columns=cols_to_drop)
isolates = bin_filt.columns

#%%
#Get isolates
isolates = data_raw.iloc[:,0]
#Format Taxa data
taxa['ID'] = taxa['ID'].str.replace('Nodule', 'N')
taxa['ID'] =  taxa['ID'].str.replace('Root', 'R')
taxa.rename(columns= {'ID':'Sample'}, inplace = True)

# # Extract day and batch Information
days1 = []
batches1 = []
for cols in data_raw.columns:
    if cols.startswith('DPI'):
        components = cols.split("_")
        batches = components[2]
        days = components[1]
        days1.append(days)
        batches1.append(batches)
           
    else:
        pass

# #%%
# #Generated better formatted columns and set batch and time points
columnHeadersClean = data_raw.columns[1:]
batches = set(batches1)
batches = list(batches)

#%%
### Create a list of DFs by batch - curently only being used for dynGenie formatting
batchDFs = []
for i in range(len(batches)):
    batchDFs.append(data_raw.filter(regex=batches[i]))

# #%% Data can be discretized and tranfromed at this point for whole sample discretization or transforamtion
combinedData = pd.concat(batchDFs, axis =1, ignore_index=False)

combinedData = combinedData.T

combinedData.columns = isolates

combinedData = combinedData.T

mergedDF1 = pd.merge(combinedData, taxa, on= 'Sample', how = 'left')

genusCountsAcrossTime1 = mergedDF1.groupby('Sample')[combinedData.columns].mean()

genusCountsAcrossTime1 += 1e-6
genusCountsAcrossTime1 = genusCountsAcrossTime1.T

genusCountsIndex = genusCountsAcrossTime1.index
genusCountsColumns = genusCountsAcrossTime1.columns


# #%% test
genusCountsAcrossTime2 = genusCountsAcrossTime1.to_numpy().flatten()
threshold = np.percentile(genusCountsAcrossTime2, 95)
# plt.hist(genusCountsAcrossTime2,bins=150)
# # Add a vertical red line at x=700
# plt.axvline(x=threshold, color='red', linestyle='--', label='x = 707.11')

# # Customize the appearance of the line (optional)
# plt.text(threshold, 10, 'x = 707.11', rotation=90, color='red')

# # Set labels and title
# plt.xlabel('ASV Count')
# plt.ylabel('Frequency')
# plt.title('Count distribution with Red Vertical Line the 95th Percentile')

# # Show the legend
# plt.legend()

# # Show the plot
# plt.show()

#%%
# Set a threshold
# threshold = threshold

# # Create a boolean mask that identifies values above or equal to the threshold
# mask = genusCountsAcrossTime2 >= threshold

# # Use the mask to index the array and keep values above or equal to the threshold
# filtered_arr = genusCountsAcrossTime2[mask]

# plt.hist(filtered_arr,bins=10)
# # Add a vertical red line at x=700
# # plt.axvline(x=threshold, color='red', linestyle='--', label='x = 707.11')

# # Customize the appearance of the line (optional)
# # plt.text(threshold, 10, 'x = 707.11', rotation=90, color='red')

# # Set labels and title
# plt.xlabel('ASV Count')
# plt.ylabel('Frequency')
# plt.title('Count distribution after removal of values beneath 95%')

# Show the legend
# plt.legend()

# # Show the plot
# plt.show()

# #%% Thresholding
# result = genusCountsAcrossTime1.applymap(lambda x: 1 if x > threshold else 0)


#%% The below code can be run to calculate the columns that contain the most variance (by isolate currently) and then creates a dataframe with them. 

# variances = genusCountsAcrossTime1.var()
# sorted_var = variances.sort_values(ascending=False)
# top_50_columns = sorted_var[:50]
# genusCountsAcrossTime1 = genusCountsAcrossTime1[top_50_columns.index]
# genusCountsIndex = genusCountsAcrossTime1.index
# genusCountsColumns = genusCountsAcrossTime1.columns

#%% Generate a random matrix of numbers # Use this if wanting to generate a random matrix of values
# genusCountsAcrossTime1 = genusCountsAcrossTime1.T
# mean_ran = genusCountsAcrossTime1.mean()
# std_dev = genusCountsAcrossTime1.std()
# min_value = genusCountsAcrossTime1.min()
# max_value = genusCountsAcrossTime1.max()

# # random_matrix = np.random.uniform(min_value, max_value, size=genusCountsAcrossTime1.shape)

# random_matrix = np.random.normal(loc=mean_ran, scale=std_dev, size=genusCountsAcrossTime1.shape)
# random_matrix = np.abs(random_matrix)
# random_df = pd.DataFrame(random_matrix, columns=genusCountsAcrossTime1.columns, index=genusCountsAcrossTime1.index)

# genusCountsAcrossTime1 = random_df


# #%% Turn rows into numpy Arrays for dynGenie
# final_array = []
# for i in batchDFs:
#     time_array = []
#     for q in i.columns:
#         i[q] = (i[q]/sum(i[q]))
#     transformed = i.T
#     for j in range(len(transformed)):
#         row_as_array = transformed.iloc[j].values
#         time_array.append(row_as_array)
#     final_array.append(time_array)

# #%% Time series Arrays for DynGenie

# timepoint1 = np.array([1,2,3,5,7,10,21,28,35,42,63])
# timepoint2 = np.array([1,2,3,5,7,10,21,28,35,42,63])
# timepoint3 = np.array([1,2,3,5,7,10,21,28,35,42,63])
# timepoint4 = np.array([1,2,3,5,7,10,21,28,35,42,63])

# final_tp = [timepoint1, timepoint2, timepoint3, timepoint4]

#%% This is where data manipulation can occur for transformations/normalizations - check thos

# TSS normalizationm - the mean after TSS is used as input for discretization
# genusCountsAcrossTime1 = genusCountsAcrossTime1
# for i in genusCountsAcrossTime1.columns:
#     genusCountsAcrossTime1[i] = (genusCountsAcrossTime1[i]/sum(genusCountsAcrossTime1[i]))*100

# genusMax = genusCountsAcrossTime1.max(axis=0)
# genusMin = genusCountsAcrossTime1.min(axis=0)
# genusMean = genusCountsAcrossTime1.mean(axis=0)
# genusMedian = genusCountsAcrossTime1.median(axis=0) # Compare to qcut binning

# genusCountsAcrossTime1 = genusCountsAcrossTime1.T
#%%CLR normalization - could use the mean of CLR for mean clr Transformation
# genusCountsAcrossTime1 = clr(genusCountsAcrossTime1)
# genusCountsAcrossTime1 = pd.DataFrame(genusCountsAcrossTime1, columns = genusCountsColumns, index= genusCountsIndex)
# genusCountsAcrossTime1 = genusCountsAcrossTime1.T

# genusMax = genusCountsAcrossTime1.max(axis=0)
# genusMin = genusCountsAcrossTime1.min(axis=0)
# genusMean = genusCountsAcrossTime1.mean(axis=0)
# genusMedian = genusCountsAcrossTime1.median(axis=0) # Compare to qcut binning


#%% Test
# genusCountsAcrossTime1 = genusCountsAcrossTime1.to_numpy().flatten()

#%% binning
# bins = pd.concat([genusMin, genusMean, genusMax], axis=1)
# bins = bins.T
# #Currently discritized by sample within batch (Could think about discritizing by batch or isolate across time within batch)
# ### Need a mean discretizer in addition to the median one provided - set numberofbins to x and measure='median' for median calc, for mean set numberbins=bins and measure='mean'
# # final_decomp1, binRanges1 = binningZeroOne(genusCountsAcrossTime1, genusCountsIndex, 2,[0,1], 'median')

# genusCountsAcrossTime1 = genusCountsAcrossTime1.T # This line is run if you are discretizing across isolates

# quantile_decomp, qCut_bin_ranges = quantile_cut(genusCountsAcrossTime1)
# subset_df = quantile_decomp.loc[:, [col for col in quantile_decomp.columns if '_bins' in col]]

# subset_df.columns = subset_df.columns.str.replace('_bins', '')

# # These lines can be run if you are discretizing aross isolates
# subset_df = subset_df.T


# #%% Remove values after discretization that don't change - if no using this ensure that isolates is checked below on all code - also check for final_decomp1
# removingNonChangers = result
# unique_values = removingNonChangers.nunique()
# constant_columns = unique_values[unique_values == 1].index.tolist()
# removingNonChangers = removingNonChangers.drop(constant_columns, axis=1)
# final_decomp2 = removingNonChangers
# isolates1 = final_decomp2.columns



#%% Removing values less than 1%
# genusCountsAcrossTime2 = genusCountsAcrossTime1[genusCountsAcrossTime1 >= 1]
# nanCount = genusCountsAcrossTime2.isnull().sum(axis=1)

# for i in genusCountsAcrossTime1.columns:
#     genusCountsAcrossTime1[i] = genusCountsAcrossTime1.loc[genusCountsAcrossTime1[i] >= 1 ]
#%%Clustering Methods - something to consider in the near future
# clustering = AgglomerativeClustering(n_clusters = None,linkage='average', compute_distances=True, distance_threshold=4)

# clustering1 = clustering.fit_predict(genusCountsAcrossTime1)
# distance = clustering.distances_

# dbscan_Clustering = DBSCAN(eps=3, min_samples=4).fit_predict(subset_df)
# unique_count = Counter(dbscan_Clustering)
# print(unique_count)

# kmeansClustering = KMeans(n_clusters=10, n_init=100).fit_predict(genusCountsAcrossTime1)

#%% Get data sorted by day and extract the days - for original data and data created based on a guassian curve of all data
# f = subset_df

# isolates1 = f.index
# f = f.T
data = bin_filt
data['Day'] = data.index.str.split('_').str[1]
data['Day'] = pd.to_numeric(data['Day'])
data['Rep'] = data.index.str.split('_').str[2]
data['Rep'] = pd.to_numeric(data['Rep'])
data = data.sort_values(by=['Day', 'Rep'])
days = sorted(list(set(data['Day'])))
reps = sorted(list(set(data['Rep'])))

#%% get the data sorted by day and rep based on multiple string values
data = bin_filt
split_index = data.index.str.split('_')
data['Day'] = split_index.str[1]
data['Day'] = pd.to_numeric(data['Day'])
data['Rep'] = split_index.str[2:].str.join('_')
# data = data.sort_values(by=['Day'])
days = sorted(list(set(data['Day'])))
reps = sorted(list(set(data['Rep'])))

#%% Get batch level data for confusion matrix comparisons

# batch1DF = final_decomp2[final_decomp2['Rep'] == 'B1']
# batch1DF = batch1DF.drop(['Rep', 'Day'], axis=1)
# batch1DF.to_csv('
# batch1DF = batch1DF.pop(batch1DF.columns[-2])

#%%
# ### Get the unique data windows and columns associated
dataWindowDay1 = []
dataWindows = []
bnstructHeaders1 = []
for i in range(len(days)-1):

    cond1 = data['Day'] >= days[i]
    cond2 = data['Day'] <= days[i+1]
    dbn_samples = data[cond1 & cond2]
    # dbn_samples = data.index[data['Day'] >= days[i]] & data.index[data['Day'] <= days[i+1]]

    datawindowDay = []
    for i in dbn_samples.index:
        dayiso = i.split('_')
        datawindowDay.append(int(dayiso[1]))

    datawindowDay = sorted(list(set(datawindowDay)))
    print(datawindowDay)
    bnstructHeaders = []
    for i in datawindowDay:
        for j in isolates:
            bnstructHeaders.append(j + '_' + str(i))

    dataWindowDay1.append(datawindowDay)
    dataWindows.append(dbn_samples)
    bnstructHeaders1.append(bnstructHeaders)
#%% Get the data that is associated with the different windowed samples
bnstructWindow = []
for i in range(len(dataWindows)):
    if dataWindows[i] is not None:
        t0_data = data.loc[dataWindows[i].index]
        bnstructWindow.append(t0_data)


#%% Bnstruct data creation across all windows in proper batch format regardless of size
timewindows = []
for i in bnstructWindow:
    li = []
    for j in reps:
        b1 = i[i['Rep'] == j]
        # print(b1)
        b1 = b1.drop(columns=['Day', 'Rep'])
        # print(b1)
        b1 = b1.to_numpy().flatten()
        # print(b1)
        li.append(b1)
        
    arry = np.array(li)
    df = pd.DataFrame(arry)
    df += 1
    timewindows.append(df)

#%%
### create data for BNStruct Across all windows in proper batch format
# timewindows = []
# for i in bnstructWindow:
#     b1 = i[i['Rep'] == 1]
#     b1 = b1.drop(columns=['Day', 'Rep'])
#     b1 = b1.to_numpy().flatten()

#     b2 = i[i['Rep'] == 2]
#     b2 = b2.drop(columns=['Day', 'Rep'])
#     b2 = b2.to_numpy().flatten()

#     b3 = i[i['Rep'] == 4]
#     b3 = b3.drop(columns=['Day', 'Rep'])
#     b3 = b3.to_numpy().flatten()

#     b4 = i[i['Rep'] == 4]
#     b4 = b4.drop(columns=['Day', 'Rep'])
#     b4 = b4.to_numpy().flatten()

#     li = [b1, b2, b3, b4]
#     arry = np.array(li)
#     df = pd.DataFrame(arry)
#     df += 1
#     timewindows.append(df)

#%% Create the bulk out files for 2tsBN
mergedData = pd.concat(timewindows)

batch1 = mergedData.loc[0]
batch2 = mergedData.loc[1]
batch3 = mergedData.loc[2]
batch4 = mergedData.loc[3]

#%% write out the bulk files
mergedData.to_csv('all_data.txt', sep=' ', header=False, index=False)
batch1.to_csv('/isolate_batch1_95.txt', sep=' ', header=False, index=False)
batch2.to_csv('/isolate_batch2_95.txt', sep=' ', header=False, index=False)
batch3.to_csv('/isolate_batch3_95.txt', sep=' ', header=False, index=False)
batch4.to_csv('/isolate_batch4_95.txt', sep=' ', header=False, index=False)

#%%
### Gets layers for data frames for BN Struct
layersForHeader = []
for i in bnstructHeaders1:
    timewindowlayer = []
    for j in i:
        day = j.split('_')
        timewindowlayer.append(day[1])
    layersForHeader.append(timewindowlayer)

layers1 = []
for i in range(len(dataWindowDay1[0])):
    for j in range(len(isolates)):
        layers1.append(i)

layers1 = [x+1 for x in layers1]
 
#%% Get headers for headers file input to BNStruct
discritization = []
state = []

for i in layersForHeader:
    discritization1 = []
    state1 = []
    for j in range(len(i)):
        discritization1.append('d')
        state1.append(2)
    discritization.append(discritization1)
    state.append(state1)

listOfHeaderFiles = []
for i in range(len(bnstructHeaders1)):
    headerFile = pd.DataFrame(columns=bnstructHeaders1[i])
    headerFile.loc[len(headerFile)] = state[i]
    headerFile.loc[len(headerFile)] = discritization[i]
    listOfHeaderFiles.append(headerFile)
# #%%
# ### Write data for BNStruct (windows of data)
# for i in range(len(timewindows)):
#     timewindows[i].to_csv('timewindow_' + str(i) + '.txt', index=False, header=False, sep = ' ')

#%%
# ### Write header file for BNStruct (windows of data)
# for i in range(len(listOfHeaderFiles)): 
#     listOfHeaderFiles[i].to_csv('r_' + str(i) + '.txt', index=False, sep = " ", quoting=csv.QUOTE_NONNUMERIC, escapechar='\\')

#%% Write single header file for bulk data (t2bn) #Issue here with quotes around the discretizing variable"
listOfHeaderFiles[0].to_csv('isolate_complete_header.txt', index=False, sep = " ", quoting=csv.QUOTE_NONNUMERIC, escapechar='\\')

#%%
### Write layer file for BNStruct (a single layer file)
with open(r'/isolate_complete_layers.txt', 'w') as fp:
    for i in layers1:
    # write each item on a new line
        fp.write("%s\n" % i)
    print('Done')

# %%
