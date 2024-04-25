## Notes
# Definitions
# CIRI2_#junction_reads = ciri2 read counts
# CLEAR_readNumber = clear read counts

## Things to do with this data
# 1 filter to runs where only deep or shallow sequencing reads are compared - no need to TPM normalize if there is no combing across data
# 2 remove reads at 1 or below
# 3 analyze data from CIRI2 or CLEAR in isolation in either deep or shallow sequencing
# 4 seperate data by round and compare as there seems to be high correlation amongts genes from certain rounds

#%% Read in your libraries
import pandas as pd
import numpy as np
import seaborn as sns
import PyWGCNA
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

#%% Read in your files circRNA Files
lotusCircRNA = pd.read_csv('data.csv', index_col=False)

#%% Create a unique column for futher processing as well as replace the comma seperating id names witha  semicolon
lotusCircRNA['Unique'] = lotusCircRNA['Sample_number'].astype(str) + '_' + lotusCircRNA['Depth']

lotusCircRNA['New_ID'] = lotusCircRNA['New_ID'].str.replace(',', ';')

#%% Create deep dataframe
deepDF = lotusCircRNA.loc[lotusCircRNA['Depth'] == 'Deep']

#%% Create a CLEAR Dataframe
clearDF = deepDF.loc[deepDF['Analysis'] == 'CLEAR']
deepDF = clearDF

#%% Create a Clear Dataframe with only roud D data
dDF = lotusCircRNA.loc[lotusCircRNA['Round'] == 'd']
deepDF = dDF

#%% Create a CIRI Dataframe
# ciriDF = lotusCircRNA.loc[lotusCircRNA['Analysis'] == 'CIRI2']
# deepDF = ciriDF

#%% Combine CLEAR and CIRI2 reads into a single column
deepDF['Total_Count'] = deepDF['CLEAR_readNumber'].add(deepDF['CIRI2_#junction_reads'], fill_value=0)

#%% Create deep pivot table based on isolated data from Ciri and Clear data
lotusCircDeep = pd.pivot_table(deepDF, values = 'Total_Count', index = 'Unique', columns='New_ID', aggfunc='sum', fill_value=0)

lotusCircDeep.rename(columns=lambda x: x.replace(',', ';'), inplace=True)

#%% Create a deep pivot table based on isolate data from Clear Data
# lotusCircDeep = pd.pivot_table(deepDF, values = 'CLEAR_readNumber', index = 'Unique', columns='New_ID', aggfunc='sum', fill_value=0)

# lotusCircDeep.rename(columns=lambda x: x.replace(',', ';'), inplace=True)

# #%% Create a deep pivot table based on isolate data from CIRI2 Data
# lotusCircDeep = pd.pivot_table(deepDF, values = 'CIRI2_#junction_reads', index = 'Unique', columns='New_ID', aggfunc='sum', fill_value=0)

# lotusCircDeep.rename(columns=lambda x: x.replace(',', ';'), inplace=True)

#%% Drop reads for genes that are 4 or less
column_sums = lotusCircDeep.sum()
columns_to_drop =  column_sums[column_sums <= 1].index
lotusCircDeep = lotusCircDeep.drop(columns=columns_to_drop)

#%% create Sample Meta for symbiosis Experiment
sample_meta = deepDF[['Unique', 'Treatment', 'Round', 'Depth']]
sample_meta = sample_meta.drop_duplicates(subset=['Unique'], keep='first')
sample_meta = sample_meta.reset_index(drop=True)
sample_meta.set_index('Unique', inplace=True)

#%% Create Gene_Meta_file for symbiosis Experiment
geneMeta = lotusCircDeep.T
colsToDrop = geneMeta.columns.to_list()

geneMeta = geneMeta.drop(columns=colsToDrop)
lotusCircRNAFiltered = lotusCircRNA.drop_duplicates(subset =['New_ID'], keep='first')

lotusCircRNAFiltered = lotusCircRNAFiltered.drop(columns=['chr', 'start', 'end', 'Round',
       'Sample_number', 'Treatment', 'Analysis', 'Tissue', 'Depth',
       'circRNA_ID', 'Repeat_overlap',
       'CIRI2_#junction_reads', 'CIRI2_SM_MS_SMS', 'CIRI2_#non_junction_reads',
       'CIRI2_junction_reads_ratio', 'CIRI2_gene_id2', 'strand',
       'CIRI2_junction_reads_ID', 'CIRI2_Product2', 'CLEAR_name',
       'CLEAR_score', 'CLEAR_exonCount', 'CLEAR_exonSizes',
       'CLEAR_exonOffsets', 'CLEAR_readNumber', 'CLEAR_isoformName',
       'CLEAR_index', 'CLEAR_flankIntron', 'CLEAR_FPBcirc', 'CLEAR_FPBlinear',
       'CLEAR_CIRCscore', 'Unique'])

mergedGeneMeta = geneMeta.merge(lotusCircRNAFiltered, how ='left', left_index=True, right_on='New_ID')
mergedGeneMeta.set_index('New_ID', inplace=True)

#%% Begin object creation for pyWGCNA
allGeneExp = PyWGCNA.WGCNA(name = 'deep', species='bacteria', geneExpPath='gene_exp.csv', outputPath='output', save=True)

#%% Perform initial analysis
allGeneExp.preprocess()
#%%
allGeneExp.findModules()
#%%
allGeneExp.updateSampleInfo(path = 'sample_meta.csv', sep = ',')

allGeneExp.setMetadataColor('Depth', {'Deep': 'deeppink', 'Shallow': 'darkviolet'})
allGeneExp.setMetadataColor('Treatment', {'C': 'thistle', 'M': 'violet', 'R': 'purple', 'A': 'magenta'})
# allGeneExp.setMetadataColor('Round', {'d': 'green', 'b': 'yellow', 'c': 'red', 'a': 'blue'})

#%%
allGeneExp.updateGeneInfo(path = 'gene_meta.csv', sep = ',')

#%%
allGeneExp.analyseWGCNA()
# allGeneExp.top_n_hub_genes(moduleName="brown", n=100)
#%%
allGeneExp.saveWGCNA()

# Initial analysis completed

# %%
