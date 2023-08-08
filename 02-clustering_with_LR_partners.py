#!/usr/bin/env python

'''
This script assesses whether by considering ligand-receptor pairs that are disrupted in diseases we have substantially more overlaps in risk across disorders. 

INPUTS: 
    - OmnipathDB ligand-receptor annotation 
    - One hot encoded matrix of genes v diseases
PROCESS:
    - Fetch all LR interactions where at least one interactor associated with disease
    - Log these interactions in a dataframe, each LR interaction v diseases 
    - Calculate adjacency matrix showing number of overlapping interactions across diseases

OUTPUTS: 
    - Heatmap of disease overlaps 
    - Adjacency matrix with couunts of overlapping LR interactions

'''

import omnipath as op
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pySankey import sankey as sn
import mellon as ml

db = op.interactions.import_intercell_network(transmitter_params = {"categories":"ligand"}, receiver_params = {"categories": "receptor"})
db = db[np.logical_not(db['genesymbol_intercell_source'].str.startswith('HLA'))]
db = db[np.logical_not(db['genesymbol_intercell_target'].str.startswith('HLA'))]
db = db[~db['genesymbol_intercell_target'].astype(str).str.startswith('COMPLEX')]
db = db[~db['genesymbol_intercell_source'].astype(str).str.startswith('COMPLEX')]
tgts_all = pd.read_csv('/Users/melis/Documents/Gene-targets/Disease_gene_associations_OpenTargets/hot-encoded-diseases.csv')

tgts_all['gene']=tgts_all['genes']
tgts = tgts_all.drop(columns=['genes']).set_index(['gene'])
tgts

source, percSrc, targets, percTgt, sum = ml.lr.find_hits(tgts, db)
tgts.drop('Empty', axis = 1, inplace = True)

#Keep LRs where at least one interactor is associated with disease
LRs = db[['genesymbol_intercell_source','genesymbol_intercell_target']][(db['genesymbol_intercell_source'].isin(list(sum.gene))) | (db['genesymbol_intercell_target'].isin(list(sum.gene)))]

unique_genes = np.unique(list(LRs['genesymbol_intercell_source'])+list(LRs['genesymbol_intercell_target']))

#Make an empty pandas dataframe to populate with all interactions where at least one interactor is associated with disease
risk_connectors = pd.DataFrame(np.zeros(shape=(len(unique_genes), len(tgts.columns))), columns=list(tgts.columns), index=unique_genes)

#Actually populate the dataframe
for n in range(0,len(tgts.columns)):
    disease_tgts = tgts[tgts.iloc[:,n]>0.0]
    disease_source, disease_percSrc, disease_targets, disease_percTgt, disease_sum = ml.lr.find_hits(disease_tgts, db)
    disease_LRs = db[['genesymbol_intercell_source','genesymbol_intercell_target']][(db['genesymbol_intercell_source'].isin(list(disease_sum.gene))) | (db['genesymbol_intercell_target'].isin(list(disease_sum.gene)))]
    disease_LRs_unique = pd.concat([disease_LRs['genesymbol_intercell_source'], disease_LRs['genesymbol_intercell_target']]).unique()
    for m in range(0,len(disease_LRs_unique)):
        print(disease_LRs_unique[m])
        print(risk_connectors.loc[disease_LRs_unique[m]])
        risk_connectors.iloc[risk_connectors.index.get_loc(disease_LRs_unique[m]), n] = 1.0

sns.clustermap(risk_connectors.T.dot(risk_connectors), vmax = 800)
plt.savefig('plots/02-clustering_with_LR_partners/clustermap_overlaps_LR_network.png', dpi = 150)
plt.show()