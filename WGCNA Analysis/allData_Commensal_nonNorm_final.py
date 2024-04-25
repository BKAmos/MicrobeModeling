## Notes
# Definitions - need to work through this document to make it suitable for the commensal data
# CIRI2_#junction_reads = ciri2 read counts
# CLEAR_readNumber = clear read counts

## Things to do with this data
# Remove leaf and root M1, P4, and V4 from the circular RNA data - this removal has been done within the Linear RNA data
# Decide on sample meta variables to consider from the Circular RNA
    # Meta - Tissue, Replicate
# Decide on Gene meta variables to consider for Circ RNA
    # CircRNA_ID, chrom, start, end, strand, gene_id, GFF3_Annotation, circRNA_type, 
# Decide on Meta variables from the linear RNA
# 
#%% Read in your libraries
import pandas as pd
import numpy as np
import seaborn as sns
import PyWGCNA
import pickle
import matplotlib.pyplot as plt

#%% Read in your files circRNA Files
lotusCircRNA = pd.read_csv('data.csv', index_col=False)


#%% Create a new column named circ len based on start and end values
lotusCircRNA['circ_len'] = lotusCircRNA['end'] - lotusCircRNA['start']

# lotusCircRNA['Unique'] = lotusCircRNA['Replicate'].astype(str) + '_' + lotusCircRNA['Tissue']
# lotusCircRNA['New_ID'] = lotusCircRNA['New_ID'].str.replace(',', ';')

#%% Drop the rows that are not of interest - these are the samples that were contaminated
# Define the specific values to drop
values_to_drop = ['LM1', 'RM1', 'LP4', 'RP4', 'LV4', 'RV4']

# Drop rows where the value in column 'Replicate' matches any of the values_to_drop
lotusCircRNA = lotusCircRNA[~lotusCircRNA['Replicate'].isin(values_to_drop)]

#%% Create a treatment columns, rename the values associated with the old nomenclature, create a round column
lotusCircRNA['Treatment'] = lotusCircRNA['Replicate'].str[1]
lotusCircRNA['New_rep'] = lotusCircRNA['Replicate'].replace({'RC1': 'RC1','RC3': 'RC2','RC4': 'RC3','RM3': 'RM1','RM4': 'RM2', 'RV1':'RMV1','RV3':'RMV2','RP1':'RMVP1','RP3':'RMVP2','LC1':'LC1','LC3': 'LC2','LC4': 'LC3','LM3': 'LM1','LM4': 'LM2', 'LV1':'LMV1','LV3':'LMV2','LP1':'LMVP1','LP3':'LMVP2'})
lotusCircRNA['Round'] = lotusCircRNA['New_rep'].str.split('([A-Za-z]+)(\d+)').str[-2]
lotusCircRNA['Treatment'] = lotusCircRNA['Treatment'].replace({'V':'MV', 'P':'MVP'})


#%% This line of code is for summing reads across dataframes for a non-normalized read count table
lotusCircRNA['Total_Count'] = lotusCircRNA['CLEAR_readNumber'].add(lotusCircRNA['CIRI2_#junction_reads'], fill_value=0)


#%% Create a Clear Dataframe with only roud X data
# lotusCircRNA = lotusCircRNA.loc[lotusCircRNA['Round'] == 'd']


#%% Create deep pivot table based on isolated data from Ciri and Clear data
lotusCircDeep = pd.pivot_table(lotusCircRNA, values = 'Total_Count', index = 'New_rep', columns='circRNA_ID', aggfunc='sum', fill_value=0)
# lotusCircDeep.rename(columns=lambda x: x.replace(',', ';'), inplace=True)


#%% Drop reads for genes based on mean read count for the frame
# all data = 1
# all data round d = 1
# all data round c = 4
# all data round b = 3
# all data round a = 1
column_sums = lotusCircDeep.sum()
columns_to_drop =  column_sums[column_sums <= 1].index
lotusCircDeep = lotusCircDeep.drop(columns=columns_to_drop)


# #%% transform the pivot table to have genes as the indexa
# lotusCircDeep = lotusCircDeep.T
# lotusCircRNA_GeneSelection = lotusCircRNA[lotusCircRNA['circRNA_ID'].isin(lotusCircDeep.index)]


# #%% The below code applies tpm normalization across two different circRNA data production files read counts. This results in a column devoid of NA and containing TPM information. Works for the symbiosis data
# lotusCircRNA_GeneSelection['CIRI2_RPK'] = (lotusCircRNA_GeneSelection['CIRI2_#junction_reads']/(lotusCircRNA_GeneSelection['circ_len']/1000))
# lotusCircRNA_GeneSelection['CLEAR_RPK'] = (lotusCircRNA_GeneSelection['CLEAR_readNumber']/(lotusCircRNA_GeneSelection['circ_len']/1000)) 
# lotusCircRNA_GeneSelection['RPK'] = lotusCircRNA_GeneSelection['CLEAR_RPK'].add(lotusCircRNA_GeneSelection['CIRI2_RPK'], fill_value=0) 
# lotusCircRNA_GeneSelection['RPK_Sum_Sample'] = lotusCircRNA_GeneSelection.groupby('New_rep')['RPK'].transform('count')
# lotusCircRNA_GeneSelection['perMillScale'] = lotusCircRNA_GeneSelection['RPK_Sum_Sample']/1000000
# lotusCircRNA_GeneSelection['TPM'] = lotusCircRNA_GeneSelection['RPK']/lotusCircRNA_GeneSelection['perMillScale']


#%% Create deep pivot table based on isolated data from Ciri and Clear data
# lotusCircDeep_2 = pd.pivot_table(lotusCircRNA_GeneSelection, values = 'TPM', index = 'New_rep', columns='circRNA_ID', aggfunc='sum', fill_value=0)
# lotusCircDeep_2.rename(columns=lambda x: x.replace(',', ';'), inplace=True)


#%% create Sample Meta for symbiosis Experiment
sample_meta = lotusCircRNA[['New_rep','Treatment', 'Tissue', 'Round']] # remove round if seperating data by round
sample_meta = sample_meta.drop_duplicates(subset=['New_rep'], keep='first')
sample_meta = sample_meta.reset_index(drop=True)
sample_meta.set_index('New_rep', inplace=True)


#%% Create Gene_Meta_file for symbiosis Experiment
geneMeta = lotusCircDeep.T
colsToDrop = geneMeta.columns.to_list()

geneMeta = geneMeta.drop(columns=colsToDrop)
lotusCircRNAFiltered = lotusCircRNA.drop_duplicates(subset =['circRNA_ID'], keep='first')

lotusCircRNAFiltered = lotusCircRNAFiltered.drop(columns=['Unnamed: 0',             'CIRI2_gene_id2', 'CIRI2_product2',
       'Experiment', 'Prep', 'Tissue', 'Replicate',
       'CIRI2_#junction_reads', 'CIRI2_SM_MS_SMS', 'CIRI2_#non_junction_reads',
       'CIRI2_junction_reads_ratio', 'CIRI2_junction_reads_ID', 'CLEAR_name',
       'CLEAR_score', 'CLEAR_exonCount', 'CLEAR_exonSizes',
       'CLEAR_exonOffsets', 'CLEAR_readNumber', 'CLEAR_isoformName',
       'CLEAR_index', 'CLEAR_flankIntron', 'CLEAR_FPBcirc', 'CLEAR_FPBlinear',
       'CLEAR_CIRCscore', 'Cluster_ID', 'circ_len', 'Treatment', 'New_rep',
       'Round', 'Total_Count'])

mergedGeneMeta = geneMeta.merge(lotusCircRNAFiltered, how ='left', left_index=True, right_on='circRNA_ID')
mergedGeneMeta.set_index('circRNA_ID', inplace=True)


#%% Begin object creation for pyWGCNA
allGeneExp = PyWGCNA.WGCNA(name = 'All data Mean Filtered', species='bacteria', geneExpPath='/gene_exp.csv', outputPath='outpath', save=True)


#%% Perform initial analysis
allGeneExp.preprocess()


#%%
allGeneExp.findModules()

#%%
allGeneExp.updateSampleInfo(path = 'sample_meta.csv', sep = ',')

allGeneExp.setMetadataColor('Tissue', {'Leaf': 'deeppink', 'Root': 'darkviolet'})
allGeneExp.setMetadataColor('Treatment', {'C': 'thistle', 'M': 'violet', 'MV': 'purple', 'MVP': 'magenta'})
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
moduleOfInterest = ['lightgrey','dimgrey','lightcoral','rosybrown','snow','gainsboro','white','whitesmoke']
# Below is the modules for round D
#['darkred','indianred','darkgrey','dimgrey','snow','brown','maroon','white','red','lightgrey','silver','firebrick','black','lightcoral']
#Below are mouldes of intererst
#Modules for all data that has been mean filtered and Normalized
#['coral', 'sienna', 'lightsalmon', 'salmon','white','rosybrown','brown', 'whitesmoke', 'darksalmon', 'lightgrey', 'snow', 'darkgrey','darkred','dimgrey','black','maroon']


for i in moduleOfInterest:
       module = pyOb.top_n_hub_genes(moduleName = i, n = 6000)
       module.to_csv('/geneList/module_' + i + '.csv')
# %%
