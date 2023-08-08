#!/usr/bin/env python

'''
This script assembles the LR network per disease and measures the distance between all networks in 
a pair-wise fashion. Outputs are summarised in a heatmap.
INPUT:
    - Omnipath database
    - One-hot encoded genes v diseases
PROCESS:
    - Fetch genes associated with all diseases
    - Filter to only include ones that act as LRs
    - Fetch from omnipath all LR interactions linked with these diseases
    - In a pairwise fashion, fetch LRs affected by each disease and measure distance between networks (DeltaCon)
    - Collate all distances in matrix
OUTPUTS: 
    - .csv with all distances
    - heatmap visualising distances between diseases
'''

import omnipath as op
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pySankey import sankey as sn
import mellon as ml
import networkx as nx
from matplotlib.pyplot import figure
import netrd
import itertools

db = op.interactions.import_intercell_network(transmitter_params = {"categories":"ligand"}, receiver_params = {"categories": "receptor"})
db = db[np.logical_not(db['genesymbol_intercell_source'].str.startswith('HLA'))]
db = db[np.logical_not(db['genesymbol_intercell_target'].str.startswith('HLA'))]
db = db[~db['genesymbol_intercell_target'].astype(str).str.startswith('COMPLEX')]
db = db[~db['genesymbol_intercell_source'].astype(str).str.startswith('COMPLEX')]


def make_disease_network(unique_genes, tgts, dis1, dis2):   
    #Make this a function in the library

    source, percSrc, targets, percTgt, sum = ml.lr.find_hits(tgts[tgts[dis1]>0], db)
    LRs = db[['genesymbol_intercell_source','genesymbol_intercell_target']][(db['genesymbol_intercell_source'].isin(list(sum.gene))) | (db['genesymbol_intercell_target'].isin(list(sum.gene)))]

    # Adjacency matrix with entire network across all diseases
    adj_LRs = pd.DataFrame(np.zeros(shape=(len(unique_genes), len(unique_genes))), index = unique_genes, columns = unique_genes)

    for n in range (0, np.shape(LRs)[0]):
        adj_LRs.iloc[adj_LRs.index.get_loc(LRs.iloc[n]['genesymbol_intercell_source']), adj_LRs.columns.get_loc(LRs.iloc[n]['genesymbol_intercell_target'])] += 1

    #Number of cell-cell interactions ligands of interest are involved in
    print(adj_LRs.sum(axis = 0).sort_values(ascending=False)[adj_LRs.sum(axis = 0).sort_values(ascending=False).index.isin(list(tgts.index))])

    #Number of cell-cell interactions receptors of interest are involved in
    adj_LRs.sum(axis = 1).sort_values(ascending=False)[adj_LRs.sum(axis = 1).sort_values(ascending=False).index.isin(list(tgts.index))]

    # Overview of LR network and calculating rank of importance with PageRank

    G1 = nx.from_pandas_adjacency(adj_LRs)

    source, percSrc, targets, percTgt, sum = ml.lr.find_hits(tgts[tgts[dis2]>0], db)
    LRs = db[['genesymbol_intercell_source','genesymbol_intercell_target']][(db['genesymbol_intercell_source'].isin(list(sum.gene))) | (db['genesymbol_intercell_target'].isin(list(sum.gene)))]

    # Adjacency matrix with entire network across all diseases
    adj_LRs = pd.DataFrame(np.zeros(shape=(len(unique_genes), len(unique_genes))), index = unique_genes, columns = unique_genes)

    for n in range (0, np.shape(LRs)[0]):
        adj_LRs.iloc[adj_LRs.index.get_loc(LRs.iloc[n]['genesymbol_intercell_source']), adj_LRs.columns.get_loc(LRs.iloc[n]['genesymbol_intercell_target'])] += 1

    #Number of cell-cell interactions ligands of interest are involved in
    print(adj_LRs.sum(axis = 0).sort_values(ascending=False)[adj_LRs.sum(axis = 0).sort_values(ascending=False).index.isin(list(tgts.index))])

    #Number of cell-cell interactions receptors of interest are involved in
    adj_LRs.sum(axis = 1).sort_values(ascending=False)[adj_LRs.sum(axis = 1).sort_values(ascending=False).index.isin(list(tgts.index))]

    # Overview of LR network and calculating rank of importance with PageRank

    G2 = nx.from_pandas_adjacency(adj_LRs)

    delta = netrd.distance.DeltaCon()
    distance = delta.dist(G1, G2)
    return G1, G2, distance

######------------------------------------
all_res = ['0.700000', '0.400000', '0.100000']
for res in range(0,len(all_res)):
    tgts_all = pd.read_csv('processed_data/hot-encoded-diseases_'+all_res[res]+'.csv')
    tgts_all['gene']=tgts_all['genes']
    tgts = tgts_all.drop(columns=['genes']).set_index(['gene'])

    source, percSrc, targets, percTgt, sum = ml.lr.find_hits(tgts, db)

    LRs = db[['genesymbol_intercell_source','genesymbol_intercell_target']][(db['genesymbol_intercell_source'].isin(list(sum.gene))) | (db['genesymbol_intercell_target'].isin(list(sum.gene)))]

    unique_genes = np.unique(list(LRs['genesymbol_intercell_source'])+list(LRs['genesymbol_intercell_target']))

    tgts = tgts.drop(columns = ['Empty'])

    dis = list(itertools.combinations(list(tgts.columns), 2))
    deltacon_results = pd.DataFrame(np.zeros(shape=(len(tgts.columns), len(tgts.columns))), index = list(tgts.columns), columns = list(tgts.columns))

    for n in range(0,len(dis)):
        G1, G2, distance = make_disease_network(unique_genes, dis[n][0], dis[n][1])

        print(dis[n][0])
        print(dis[n][1])

        deltacon_results.loc[dis[n][0]][dis[n][1]] = distance
        deltacon_results.loc[dis[n][1]][dis[n][0]] = distance
        deltacon_results.to_csv('processed_data/03-LR_network_visualisation/deltacon_results_crossdisease'+all_res[res]+'.csv')

    sns.clustermap(deltacon_results)
    plt.savefig(('plots/03-LR_network_visualisation/deltacon_results_crossdisease_'+all_res[res]+'thr.csv'))
    plt.show()