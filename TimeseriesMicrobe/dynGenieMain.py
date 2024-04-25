#%%
from dynGENIE3 import *
import _pickle
import numpy as np
import pandas as pd


#%%
### import the data
data_raw = pd.read_csv('FilteredData.csv')
taxa = pd.read_csv('/taxonomy.csv', sep = '\t')
meta = pd.read_csv('ES_Meta.csv')

#%%
#Gather Isolate Data
isolates = data_raw.iloc[:,0]

#Format Taxa data
taxa['ID'] = taxa['ID'].str.replace('Nodule', 'N')
taxa['ID'] =  taxa['ID'].str.replace('Root', 'R')
taxa.rename(columns= {'ID':'Sample'}, inplace = True)

# Extract day and batch Information
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

#%%
#Generated better formatted columns and set batch and time points
columnHeadersClean = data_raw.columns[1:]
batches = set(batches1)
batches = list(batches)

#%%
### Create a list of DFs by batch - curently not using
batchDFs = []
for i in range(len(batches)):
    batchDFs.append(data_raw.filter(regex=batches[i]))

### Data can be discretized and tranfromed at this point for whole sample discretization or transforamtion
combinedData = pd.concat(batchDFs, axis =1, ignore_index=False)

combinedData = combinedData.T

combinedData.columns = isolates

combinedData = combinedData.T

mergedDF1 = pd.merge(combinedData, taxa, on= 'Sample', how = 'left')

genusCountsAcrossTime1 = mergedDF1.groupby(['Sample']).mean()

genusCountsAcrossTime1 += 1e-6
genusCountsAcrossTime1 = genusCountsAcrossTime1.T

genusCountsIndex = genusCountsAcrossTime1.index
genusCountsColumns = genusCountsAcrossTime1.columns

#%% Turn rows into numpy Arrays for dynGenie
final_array = []
for i in batchDFs:
    # time_array = []
    for q in i.columns:
        i[q] = (i[q]/sum(i[q]))
    transformed = i.T
    print(transformed)
    single_array = transformed.to_numpy()
    final_array.append(single_array)
    # for index,row in transformed.iterrows():
    #     row_as_array = np.array(row)
    #     final_array.append(row_as_array)
    # final_array.append(time_array)

#%% Time series Arrays for DynGenie

timepoint1 = np.array([1,2,3,5,7,10,21,28,35,42,63])
timepoint2 = np.array([1,2,3,5,7,10,21,28,35,42,63])
timepoint3 = np.array([1,2,3,5,7,10,21,28,35,42,63])
timepoint4 = np.array([1,2,3,5,7,10,21,28,35,42,63])

final_tp = [timepoint1, timepoint2, timepoint3, timepoint4]

#%%
isolates = isolates.to_list()
# %% Below is the test code that comes with the tutorial
# f = open('TS_data.pkl','rb')
# (TS_data, time_points, decay_rates, gene_names) = _pickle.load(f)
# f.close()
# # %%
# (VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(TS_data, time_points, compute_quality_scores=True, save_models=True)
# %% Test Data
(VIM, alphas, prediction_score, stability_score, treeEstimators) = dynGENIE3(TS_data = final_array, time_points= final_tp, gene_names=isolates, compute_quality_scores=True, save_models=True, nthreads=1)

#%%
interactions = get_link_list(VIM, gene_names=isolates)

# %%
help(dynGENIE3)
# %%
