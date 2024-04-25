#%%
# Linear RNA analysis to Do
# remove "_S398_L003.x.clean.bam" from column names - the _S398 is unique to each column
# Drop "Product" column
# Calculate mean and remove genes with counts beneath that are - gene exp 
# Gene_meta - get from geneID and Product column
# Sample Meta - generate based on transposed gene exp


#%% Load the librarties to use
import pandas as pd
import PyWGCNA


#%% Read in the linear RNA files
lotusLinRNA = pd.read_csv('data.csv', index_col=False)

#%% Define a function to determine the tissue based on the 'Sample' column value
def get_tissue(sample):
    if 'L' in sample:
        return 'Leaf'
    elif 'R' in sample:
        return 'Root'
    else:
        return None  # Or any other value if the condition is not met
    

#%% Rename columns by removing strings after the second '_'
new_columns = []
for column in lotusLinRNA.columns:
    if '_' in column:
        parts = column.split('_')
        new_column = '_'.join(parts[:2])
        new_column = new_column.replace('_', '')
        new_columns.append(new_column)
    else:
        new_columns.append(column)

lotusLinRNA.columns = new_columns

#%% Isolate gene info
geneInfo = lotusLinRNA[['geneId', 'Product']]


# %% Drop geneID and Product, transpose, create a new column with correct IDs for sample and sample Meta
lotusLinRNASamp = lotusLinRNA.drop(columns=['Product'])
lotusLinRNASamp.set_index('geneId', inplace=True)
lotusLinRNASamp = lotusLinRNASamp.rename_axis(None)
lotusLinRNASamp = lotusLinRNASamp.T


#%% Below processign can occur after filtering
lotusLinRNASamp.reset_index(inplace=True)
lotusLinRNASamp = lotusLinRNASamp.rename(columns={'index':'Sample'})
lotusLinRNASamp['Treatment'] = lotusLinRNASamp['Sample'].str[1]
lotusLinRNASamp['New_rep'] = lotusLinRNASamp['Sample'].replace({'Rc1': 'Rc1','Rc3': 'Rc2','Rc4': 'Rc3','Rm3': 'Rm1','Rm4': 'Rm2', 'Rv1':'Rmv1','Rv3':'Rmv2','Rp1':'Rmvp1','Rp3':'Rmvp2','Lc1':'Lc1','Lc3': 'Lc2','Lc4': 'Lc3','Lm3': 'Lm1','Lm4': 'Lm2', 'Lv1':'Lmv1','Lv3':'Lmv2','Lp1':'Lmvp1','Lp3':'Lmvp2'})
lotusLinRNASamp['Round'] = lotusLinRNASamp['New_rep'].str.split('([A-Za-z]+)(\d+)').str[-2]
lotusLinRNASamp['Treatment'] = lotusLinRNASamp['Treatment'].replace({'v':'mv', 'p':'mvp'})
lotusLinRNASamp['Tissue'] = lotusLinRNASamp['Sample'].apply(get_tissue)


#%% create Sample Meta for commensal Experiment
sample_meta = lotusLinRNASamp[['New_rep','Treatment', 'Tissue', 'Round']] # remove round if seperating data by round
sample_meta = sample_meta.drop_duplicates(subset=['New_rep'], keep='first')
sample_meta = sample_meta.reset_index(drop=True)
sample_meta.set_index('New_rep', inplace=True)
sample_meta = sample_meta.rename_axis(None)


#%%
lotusLinRNASamp.set_index('New_rep', inplace=True)
lotusLinRNASamp= lotusLinRNASamp.rename_axis(None)
lotusLinRNASamp = lotusLinRNASamp.drop(columns=['Sample','Treatment','Round','Tissue'])

#%% Drop reads for genes based on mean read count for the frame
# all data = 1463
column_sums = lotusLinRNASamp.sum()
columns_to_drop =  column_sums[column_sums <= 1463].index
lotusLinDeep = lotusLinRNASamp.drop(columns=columns_to_drop)


#%% Create Gene_Meta_file for symbiosis Experiment
genesOfInterest = lotusLinDeep.T
genesOfInterest = genesOfInterest.drop(genesOfInterest.columns, axis=1)
gene_meta = genesOfInterest.merge(geneInfo, how ='left', left_index=True, right_on='geneId')
gene_meta.set_index('geneId', inplace=True)
gene_meta = gene_meta.rename_axis(None)


# %% Begin object creation for pyWGCNA
allGeneExp = PyWGCNA.WGCNA(name = 'All linear data Mean Filtered', species='plant', geneExpPath='/gene_exp.csv', outputPath='outpath', save=True)


#%% Perform initial analysis
allGeneExp.preprocess()


#%%
allGeneExp.findModules()

#%%
allGeneExp.updateSampleInfo(path = 'sample_meta.csv', sep = ',')

allGeneExp.setMetadataColor('Tissue', {'Leaf': 'deeppink', 'Root': 'darkviolet'})
allGeneExp.setMetadataColor('Treatment', {'c': 'thistle', 'm': 'violet', 'mv': 'purple', 'mvp': 'magenta'})
allGeneExp.setMetadataColor('Round', {'1': 'green', '2': 'yellow', '3': 'red'})


#%%
allGeneExp.updateGeneInfo(path = 'gene_meta.csv', sep = ',')


#%%
allGeneExp.analyseWGCNA()
# allGeneExp.top_n_hub_genes(moduleName="brown", n=100)
#%%
allGeneExp.saveWGCNA()


# %%
pyOb = PyWGCNA.readWGCNA('data.p')
# %%
moduleOfInterest = ['dimgrey','snow']
# Below is the modules for round D
#['darkred','indianred','darkgrey','dimgrey','snow','brown','maroon','white','red','lightgrey','silver','firebrick','black','lightcoral']
#Below are mouldes of intererst
#Modules for all data that has been mean filtered and Normalized
#['coral', 'sienna', 'lightsalmon', 'salmon','white','rosybrown','brown', 'whitesmoke', 'darksalmon', 'lightgrey', 'snow', 'darkgrey','darkred','dimgrey','black','maroon']


for i in moduleOfInterest:
       module = pyOb.top_n_hub_genes(moduleName = i, n = 15000)
       module.to_csv('/geneList/module_' + i + '.csv')
# %%
