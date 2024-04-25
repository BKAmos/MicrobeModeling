#%%
import pandas as pd
from sklearn.metrics import r2_score, mean_absolute_error
import seaborn as sns

#%%
allTime = []
timeWindows = []
isolatesAcrossTime = []

#%% load the data of interest
for i in range(1,5):
    original_data = pd.read_csv('fil' + str(i) + '_g.tsv', sep='\t')

    predicted_data = pd.read_csv('file/B' + str(i) + '_g_1_predictions_mb.csv', header=None)


#isolate certain components - subject and isolate names specially
    df = original_data.drop_duplicates(subset='Subject', keep='first')
    df_index = df.index
    original_data = original_data.drop(df_index)
    original_data.reset_index(inplace=True, drop=True)
    subject = original_data['Subject']
    isolates = original_data.columns[1:]


#compare data across all time windows and batches
    original_data_for_comp = original_data.drop(columns=['Subject'])
    r2 = r2_score(original_data_for_comp, predicted_data, multioutput='variance_weighted')
    mae = mean_absolute_error(original_data_for_comp, predicted_data)
    allTime.append({'Metric': 'r2', 'Score': r2, 'Batch': i})
    allTime.append({'Metric': 'mae', 'Score': mae, 'Batch': i})
    # allTime.append({'Batch': i})


#row wise error rates (accuracy across time windows)
    for j in range(len(original_data_for_comp.index)):
        r2_time = r2_score(original_data_for_comp.loc[j], predicted_data.loc[j])
        mae_time = mean_absolute_error(original_data_for_comp.loc[j], predicted_data.loc[j])
        timeWindows.append({'Metric': 'r2', 'Score': r2_time, 'Batch': i, 'TimeWindow': j})
        timeWindows.append({'Metric': 'mae', 'Score': mae_time, 'Batch': i, 'TimeWindow': j})


#isolate wise error and accuracy
    for j in range(len(original_data_for_comp.columns)):
        r2_isolate = r2_score(original_data_for_comp.iloc[:,j], predicted_data.iloc[:,j])
        mae_isolate = mean_absolute_error(original_data_for_comp.iloc[:,j], predicted_data.iloc[:,j])
        isolatesAcrossTime.append({'Metric': 'r2', 'Score': r2_isolate, 'Batch': i, 'Isolate': original_data_for_comp.columns[j]})
        isolatesAcrossTime.append({'Metric': 'mae', 'Score': mae_isolate, 'Batch': i, 'Isolate': original_data_for_comp.columns[j]})

# %% Create dataframes

allTime1 = pd.DataFrame(allTime)
timeWindows1 = pd.DataFrame(timeWindows)
isolatesAcrossTime1 = pd.DataFrame(isolatesAcrossTime)

# %%

timeWindows1['TimeWindow'] = (timeWindows1['TimeWindow'] % 10) * (10 / (10 - 1))
timeWindows1['TimeWindow'] = timeWindows1['TimeWindow'].astype(int)
time_r2 = timeWindows1[timeWindows1['Metric']=='r2']
time_mae = timeWindows1[timeWindows1['Metric']=='mae']


#%%
times = set(time_r2['TimeWindow'].values)
for i in times:
    values = time_r2[time_r2['TimeWindow']==i]
    mean = values['Score'].mean()
    print(mean)

# timeWindows1['TimeWindow'], timeWindows1['ori'] = pd.factorize(timeWindows1['TimeWindow'])

# %%
g = sns.FacetGrid(data=timeWindows1, col='Metric')
g.map_dataframe(sns.boxplot, x='TimeWindow', y='Score')
# %%
# isolate_mapping = sns.FacetGrid(data=isolatesAcrossTime1, col='Batch', row='Metric')
# isolate_mapping.map_dataframe(sns.boxplot, x='Isolate', y='Score')
isolate_r2 = isolatesAcrossTime1[isolatesAcrossTime1['Metric']=='r2']
isolate_mae = isolatesAcrossTime1[isolatesAcrossTime1['Metric']=='mae']
# %%
