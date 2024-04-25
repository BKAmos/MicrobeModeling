#%% import libraries
import pandas as pd
from sklearn.impute import KNNImputer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

#%%
pheno = pd.read_excel('pheno.xlsx', sheet_name='Data', index_col=False)
deg = pd.read_excel('deg.xlsx')
lipid = pd.read_csv('lipid.txt', sep='\t')
lipid_names = pd.read_csv('lipid_data.tsv', sep='\t')
metabo = pd.read_csv('secondary.txt', sep = '\t')
metabo_names = pd.read_csv('metabolites_tabs.txt', sep='\t')
asv = pd.read_csv('VST_transformed_ASVs.csv')

#%% create objects to use later
imputer = KNNImputer(n_neighbors=8, weights='distance')
scaler = MinMaxScaler()

# %% Organize phenotype Data
# Drop date, phys2D, and encode the fungal variables as catagorical
pheno = pheno.drop(columns=['Date', 'PhysToD_HHMMSS'])
fungal_map = {0 : 'A',  3 : 'B',  5 : 'C', 15 : 'D', 25 : 'E', 32 : 'F', 37 : 'G', 52 : 'H', 62 : 'I'}
pheno['Fungus'] = pheno['Fungus'].map(fungal_map)

#%% One-hot-encode specific columns from the phenotypic variables
columns_to_encode = ['Block', 'Plot', 'UniqPlot', 'Subplot']
pheno_encoded = pd.get_dummies(pheno, columns=columns_to_encode)
pheno_encoded.replace({True: 1, False: 0}, inplace=True)


# %% Organize the deg Information
deg = deg.drop(columns=['Unnamed: 0'])
deg.columns = [str(col).rstrip('ABC') if str(col).endswith(('A', 'B','C')) else col for col in deg.columns]
deg.set_index('DEGenes', inplace=True)
column_mapping = {column: int(column) for column in deg.columns}
deg.rename(columns=column_mapping, inplace=True)
deg = deg.T
deg.index.name = 'UniqueNum'



# %% Organize the lipid data
# Need to pull in the names of the lipids
lipid.columns = [col.replace('Y-', '') for col in lipid.columns]
lipid.columns = [col.replace('_L_pos.mzdata', '') for col in lipid.columns]
lipid_names_1 = lipid_names['Lipid']
lipid = pd.concat([lipid_names_1, lipid], axis = 1)
lipid.set_index('Lipid', inplace=True)
column_mapping = {column: int(column) for column in lipid.columns}
lipid.rename(columns=column_mapping, inplace=True)
lipid = lipid.T
lipid.index.name = 'UniqueNum'
lipid_index = lipid.index

#%% Impute NaN Lipid data
lipid_imputed = imputer.fit_transform(lipid)
lipid_imputed = pd.DataFrame(data = lipid_imputed, index = lipid_index, columns = lipid.columns)


# %% Organize the metabolite data
# Need to pull in the names of the metabolites
metabo.columns = [col.replace('Y-', '') for col in metabo.columns]
metabo.columns = [col.replace('_S_neg_.mzdata', '') for col in metabo.columns]
metabo.columns = [col.replace('_S_neg_redo_.mzdata', '') for col in metabo.columns]
metabo_names_1 = metabo_names['Metabolite']
metabo =  pd.concat([metabo_names_1, metabo], axis = 1)
metabo.set_index('Metabolite', inplace=True)
column_mapping = {column: int(column) for column in metabo.columns}
metabo.rename(columns=column_mapping, inplace=True)
metabo = metabo.T
metabo.index.name = 'UniqueNum'
metabo_index = metabo.index

#%% Impute NaN metabolite data
metabo_imputed = imputer.fit_transform(metabo)
metabo_imputed = pd.DataFrame(data = metabo_imputed, index = metabo_index, columns = metabo.columns)


# %% Organize the ASV Data
asv.columns = [col.split('-', 1)[0] for col in asv.columns]
asv.set_index('Unnamed: 0', inplace=True)
column_mapping = {column: int(column) for column in asv.columns}
asv.rename(columns=column_mapping, inplace=True)
asv = asv.T
asv.index.name = 'UniqueNum'


# %% Merge data on the phenotype values and fill in the missing values with the KNN imputer based on certain neighbors
merge_df = pheno_encoded.merge(deg, left_on='UniqueNum', right_index=True, how='left')
merge_df = merge_df.merge(lipid_imputed, left_on='UniqueNum', right_index=True, how='left')
merge_df = merge_df.merge(metabo_imputed, left_on='UniqueNum', right_index=True, how='left')
merge_df = merge_df.merge(asv, left_on='UniqueNum', right_index=True, how='left')


#%% use KNN imputer but drop the string columns
cols_to_attach = merge_df[['Water','Fungus']].copy()
merge_df = merge_df.drop(columns=['Water','Fungus'])
merge_df.columns = merge_df.columns.astype(str)
merge_df_imputed = imputer.fit_transform(merge_df)
merge_df_imputed = pd.DataFrame(data = merge_df_imputed, index = merge_df.index, columns = merge_df.columns)
merge_df_imputed = pd.concat([merge_df_imputed, cols_to_attach], axis=1)


#%% Split the data into training and testing splits based on 80/20 splits based on fungal innoculant
# Get unique categories in the 'Fungus' column
categories = merge_df_imputed['Fungus'].unique()
watering_levels = merge_df_imputed['Water'].unique()

# Initialize empty DataFrames to store train and test data
train_data = pd.DataFrame()
test_data = pd.DataFrame()

for category in categories:
    for watering_level in watering_levels:
        # Subset data for the current category and watering level
        subset_data = merge_df_imputed[(merge_df_imputed['Fungus'] == category) & (merge_df_imputed['Water'] == watering_level)]
        
        # Split the data into train and test sets
        train_subset, test_subset = train_test_split(subset_data, test_size=0.2, random_state=42)
        
        # Append the splits to the train and test DataFrames
        train_data = pd.concat([train_data, train_subset])
        test_data = pd.concat([test_data, test_subset])


#%% Remove correlated features
# cols_to_attach = train_data[['UniqueNum','Water','Fungus']].copy()
# train_data = train_data.drop(columns=['UniqueNum','Water','Fungus'])

# threshold = 0.7  # Adjust as needed

# # Calculate the correlation matrix
# correlation_matrix = train_data.corr().abs()

# # Create a mask to identify highly correlated features
# mask = (correlation_matrix >= threshold) & (correlation_matrix < 1.0)

# # Iterate through the columns and drop one of the highly correlated columns
# columns_to_drop = set()
# for col in mask.columns:
#     correlated_cols = mask.index[mask[col]].tolist()
#     for correlated_col in correlated_cols:
#         if correlated_col not in columns_to_drop:
#             columns_to_drop.add(correlated_col)

# # Drop the highly correlated columns from the DataFrame
# train_data_low_correlation = train_data.drop(columns=columns_to_drop)

# %% Min max scale the data  
columns_to_exclude = ['UniqueNum','Water','Fungus']
# cols_to_attach = train_data[['UniqueNum','Water','Fungus']].copy()
columns_to_scale = train_data.columns.difference(columns_to_exclude)
df_scaled_train = train_data.copy()
df_scaled_train[columns_to_scale] = scaler.fit_transform(train_data[columns_to_scale])
# df_scaled_train = pd.concat([df_scaled_train, cols_to_attach], axis=1)

#%%
df_scaled_test = test_data.copy()
df_scaled_test[columns_to_scale] = scaler.fit_transform(test_data[columns_to_scale])


#%% Drop low variance columns and highly correlated columns
cols_to_attach = df_scaled_train[['Water','Fungus', 'UniqueNum']].copy()
df_scaled_train = df_scaled_train.drop(columns=['Water','Fungus', 'UniqueNum'])
#%%
var_thresh = .04
variances = df_scaled_train.var()
low_variance_columns = variances[variances < var_thresh].index
train_data_var_filtered = df_scaled_train.drop(columns=low_variance_columns)
train_data_var_filtered = pd.concat([train_data_var_filtered, cols_to_attach], axis=1)
# %% Have test dataset match training
columns_to_drop = set(df_scaled_test.columns) - set(train_data_var_filtered.columns)
df_scaled_test.drop(columns=columns_to_drop, inplace=True)

# %%
df_scaled_test.set_index('UniqueNum', inplace=True)
train_data_var_filtered.set_index('UniqueNum', inplace=True)

# %%
df_scaled_test.to_csv('varianceFiltered_MinMax_KNN_test.csv', index=True)
train_data_var_filtered.to_csv('varianceFiltered_MinMax_KNN_train.csv', index=True)
# %%
