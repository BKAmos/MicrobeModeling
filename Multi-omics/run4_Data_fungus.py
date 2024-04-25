#%% import libraries
import pandas as pd
from sklearn.model_selection import train_test_split

#%%
pheno = pd.read_excel('data.xlsx', sheet_name='Data', index_col=False)
train_run_4 = pd.read_csv('train.csv')
test_run_4 = pd.read_csv('test.csv')

# %% Drop phenotypic columns aside from the fungal treatements and convert fungal treatments into strings
pheno = pheno.drop(columns=['Block', 'Plot', 'UniqPlot', 'Water', 'Subplot',
       'Date', 'PhysToD_HHMMSS', 'Photo', 'Cond', 'Trmmol', 'WUE', 'VpdL',
       'Tair', 'Tleaf', 'SoilMoist', 'Biomass', 'Height', 'LAI'])

fungal_map = {0 : 'A',  3 : 'B',  5 : 'C', 15 : 'D', 25 : 'E', 32 : 'F', 37 : 'G', 52 : 'H', 62 : 'I'}
pheno['Fungus'] = pheno['Fungus'].map(fungal_map)
pheno.rename(columns={'UniqueNum' : 'Sample'}, inplace=True)

#%% Drop water column from training and testing data and merge based on the fungal data - This train test split is not done by fungus and is therefore a form a data leakage
# train_run_4 = train_run_4.drop(columns=['Water'])
# test_run_4 = test_run_4.drop(columns=['Water'])

train_run_fungal = pd.merge(train_run_4, pheno, on = 'Sample', how = 'left')
test_run_fungal = pd.merge(test_run_4, pheno, on = 'Sample', how = 'left')
# %% Write data leakage train/test split to files

train_run_fungal.to_csv('train_data_leakage_fungus_run5.csv', index=False)
test_run_fungal.to_csv('test_data_leakage_fungus_run5.csv', index=False)


# %% Remove data leakage by properly splitting by fungus. 80% of fungal treatement in training and 20% in testing
combined_frame = pd.concat([train_run_fungal, test_run_fungal])
combined_frame = combined_frame.reset_index(drop=True)

# %% Develope the train test splits
# Get unique categories in the 'Fungus' column
categories = combined_frame['Fungus'].unique()
watering_levels = combined_frame['Water'].unique()

# Initialize empty DataFrames to store train and test data
train_data = pd.DataFrame()
test_data = pd.DataFrame()


# # Iterate over each category and split the data
# for category in categories:
#     # Subset data for the current category
#     category_data = combined_frame[combined_frame['Fungus'] == category]
    
#     # Split the data into train and test sets
#     train_cat, test_cat = train_test_split(category_data, test_size=0.1, random_state=42)
    
#     # Append the splits to the train and test DataFrames
#     train_data = pd.concat([train_data, train_cat])
#     test_data = pd.concat([test_data, test_cat])

    # Iterate over each category and watering combination and split the data
for category in categories:
    for watering_level in watering_levels:
        # Subset data for the current category and watering level
        subset_data = combined_frame[(combined_frame['Fungus'] == category) & (combined_frame['Water'] == watering_level)]
        
        # Split the data into train and test sets
        train_subset, test_subset = train_test_split(subset_data, test_size=0.2, random_state=42)
        
        # Append the splits to the train and test DataFrames
        train_data = pd.concat([train_data, train_subset])
        test_data = pd.concat([test_data, test_subset])


#%% Write the train and test DataFrames to files
train_data.to_csv('train_data_water_fungus_split.csv', index=False)
test_data.to_csv('test_data_water_fungus_split.csv', index=False)


# %%
