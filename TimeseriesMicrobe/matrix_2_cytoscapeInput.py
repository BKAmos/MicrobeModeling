#%%
import pandas as pd

#%%
matrix = pd.read_csv('/matrix_4.txt', index_col=0, sep='\t')

#%% Get a list that contains the ending node in the i element and the starting nodes in the i+1 element
connections = []

for i in matrix.columns:
    
    j = matrix.loc[matrix[i] == 1].index.to_list()

    if len(j) >= 1:
        print(i)
        print(j)
        connections.append(i)
        connections.append(j)
        
# %% Get the starting and ending node from the matrix
startingNode = []
endingNode = []

for i in range(0,len(connections),2):

    for j in connections[i+1]:
        print(j)
        print(connections[i+1])
        # print(j)
        startingNode.append(j)
        endingNode.append(connections[i])
# %% Create a file for output to cytoscape that represents the whole DBN ###

data_tuples = list(zip(startingNode,endingNode))
to_cytoscape = pd.DataFrame(data_tuples, columns=['startingNode', 'endingNode'])

to_cytoscape.to_csv('/matrix_4_cyto.csv', index=False, sep=',')

# %% Generate two time step representations of the network
to_cytoscape['startingNode'] = to_cytoscape['startingNode'].astype(str) + '-'
to_cytoscape['endingNode'] = to_cytoscape['endingNode'].astype(str) + '-'
times = []

for i in matrix.columns:
    
    element = i.split('_')
    element = str(element[1]) + '-'
    
    if element not in times:
        times.append(element)

# %% Generate the time series DFs for analysis within Cytoscape
timeslots = []
for i in range(0,len(times)-1):
    df = to_cytoscape[((to_cytoscape['startingNode'].str.contains(times[i])) & (to_cytoscape['endingNode'].str.contains(times[i]))) | ((to_cytoscape['startingNode'].str.contains(times[i])) & (to_cytoscape['endingNode'].str.contains(times[i+1]))) | ((to_cytoscape['startingNode'].str.contains(times[i+1])) & (to_cytoscape['endingNode'].str.contains(times[i+1])))]
    timeslots.append(df)
    print(df)
    print

# %% Remove the unique idenfier so that 1 and 10 can be distinguished by .str.contains(xxxx)
for i in timeslots:
    i['startingNode'] = i['startingNode'].str.rstrip('-')
    i['endingNode'] = i['endingNode'].str.rstrip('-')

# %% Write the desired time series graphs to file
for i in range(len(timeslots)):
    timeslots[i].to_csv('/time' + str(i) +'.csv', index = False)
# %%
