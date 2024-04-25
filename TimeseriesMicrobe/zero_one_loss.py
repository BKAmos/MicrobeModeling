#%% Confusion Matrix Testing in python
import pandas as pd
import numpy as np
from sklearn.metrics import ConfusionMatrixDisplay, balanced_accuracy_score, precision_recall_fscore_support, hamming_loss, zero_one_loss
import matplotlib.pyplot as plt
import seaborn as sns

#%% Import the actual data vectors. y_true = x[n-1] xor x[n] and y_pred = y[n] xor x[n]

data_B1 = pd.read_csv('/isolate_batch1_quantile.txt', index_col=False, header=None, sep=' ')

data_B2 = pd.read_csv('isolate_batch2_quantile.txt', index_col=False, header=None, sep=' ')

data_B3 = pd.read_csv('//isolate_batch3_quantile.txt', index_col=False, header=None, sep=' ')

data_B4 = pd.read_csv('/isolate_batch4_quantile.txt', index_col=False, header=None, sep=' ')

#%% import the predicted vectors
predicted_B1 = pd.read_csv('/batch1_inf_std_noObs.txt', index_col=False, header=None)

predicted_B2 = pd.read_csv('/batch2_inf_std_noObs.txt', index_col=False, header=None)

predicted_B3 = pd.read_csv('/batch3_inf_std_noObs.txt', index_col=False, header=None)

predicted_B4 = pd.read_csv('/batch4_inf_std_noObs.txt', index_col=False, header=None)


#%% For each batch of original data select the i and i+1 column - this is the isolate level information across time, get xor data and that is the ground truth for assessment


#%% Pick rows of interest for the past data, and current data from f1 to represent x(n-1) and x(n)
def batch_iso_extraction(batch):

    batch_files = [batch]
    ori_data = []

    for i in batch_files:
        batch_data = []
        for j in range(len(i.columns)):
            batch_data.append(i.iloc[:,j])
        ori_data.append(batch_data)

    isolate_level_true_change = []
    for i in ori_data:
        for j in range(0,int(len(i)/2)):
            y_true0 = np.bitwise_xor(i[j],i[j+17])
            isolate_level_true_change.append(y_true0)
    
    return isolate_level_true_change


#%%
def class_loss_assessment(iso_level_true, iso_level_pred):
    
    class_loss_data = []
    
    for i in range(0,len(iso_level_true)):
        
        for j in range(0,len(iso_level_pred), 17):
            
            class_loss = zero_one_loss(iso_level_true[i], iso_level_pred[j+i])
            class_loss_data.append({"Isolate": i, "Classification_Loss": class_loss})

    class_loss_df = pd.DataFrame(class_loss_data)

    return class_loss_df

#%% get the predicted vector

predictions_B1 = []

# predicted_B1 = pd.concat([predicted_B1, predicted_B2, predicted_B3, predicted_B4])

for i in range(0,len(predicted_B1),34):
    
    # 46 for isoalte level random gaussian noise per batch
    # 44 for isolate level rando gaussian noise
    # 34 for isolate level 95%
    # 30 for family level non reduced
    # 10 for class level non reduced
    # 16 for order level non reduced data
    # 166 for unreduced data at isolate level
    # 22 for TSS Median Reduced isolate level
    # 126 for CLR Median Reduced isolate level
    # 84 for TSS Mean Reduced isolate level

    vector = predicted_B1.iloc[i:i+34]
    vector = vector.reset_index(drop=True)
    vector = vector.squeeze()
    predictions_B1.append(vector)


# %% subset the vector into smaller chunks and create a dataframe for each time window
dataframes = []
for i in range(0,len(predictions_B1),10):
    subset = predictions_B1[i:i+10]
    df = pd.DataFrame(subset)
    dataframes.append(df)

isolate_level_info = []
for i in dataframes:
    batch_data = []
    for j in range(len(i.columns)):
        batch_data.append(i.iloc[:,j])
    isolate_level_info.append(batch_data)


isolate_level_pred_change = []
for i in isolate_level_info:
    for j in range(0,int(len(i)/2)):
        y_true0 = np.bitwise_xor(i[j],i[j+17])
        isolate_level_pred_change.append(y_true0)


# %% MAIN
b1_true_change = batch_iso_extraction(data_B1)
b2_true_change = batch_iso_extraction(data_B2)
b3_true_change = batch_iso_extraction(data_B3)
b4_true_change = batch_iso_extraction(data_B4)

b1_class_loss = class_loss_assessment(b1_true_change, isolate_level_pred_change)
b2_class_loss = class_loss_assessment(b2_true_change, isolate_level_pred_change)
b3_class_loss = class_loss_assessment(b3_true_change, isolate_level_pred_change)
b4_class_loss = class_loss_assessment(b4_true_change, isolate_level_pred_change)
# %%
combined_dfs = pd.concat([b1_class_loss, b2_class_loss, b3_class_loss, b4_class_loss], axis=0)

# %%
