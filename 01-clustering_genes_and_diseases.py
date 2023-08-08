#!/usr/bin/env python

'''
The purpose of this script is to investigate direct overlaps in genetic risk across diseases. This is done by: 
(1) Hierarchical clustering of genes v diseases 
(2) Generating a counts matrix with number of directly overlapping genes between disorders

INPUTS: 
    - .csv file with the list of all genes and all diseases, where each element is the risk score of that gene in relation to that disease
    - .csv with hot-one-encoded genes per disease
PROCESS:
    - Performs hierarchical clustering on both axes, plots clustermap
    - Calculate counts of total number of overlapping genes across diseases 
OUTPUTS: 
    - Clustermap of genes with risk scores in relation to diseases
    - Clustermap of counts of overlapping genes between diseases
    - .csv file with adjacency matrix, where total counts of genes shared across disorders are summarised

'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

###----------------------Set some plotting parameters
sns.set(font_scale=1.1) 

#Set colours for disease type color labels
lut = {'AD' : 'purple',
 'ALS' : 'purple',
 'AnorexiaNervosa' : 'teal',
 'BipolarDisorder' : 'teal',
 'BrainAneurysm' : 'orange',
 'EssentialTremor' : 'purple',
 'FrontotemporalDementia' : 'purple',
 'IntracranialHemorrhage' : 'orange',
 'LBD' : 'purple',
 'MajorDepressiveDisorder' : 'teal',
 'MigraineDisorder' : 'orange',
 'MigraineWithAura' : 'orange',
 'MS' : 'orange',
 'NarcolepsyCataplexy' : 'orange',
 'Narcolepsy' : 'orange',
 'NeuroticDisorder' : 'teal',
 'OCD' : 'teal',
 'PartialEpilepsy' : 'orange',
 'PD' : 'purple',
 'ProgressiveSupranuclearPalsy' : 'purple',
 'RestlessLeg' : 'orange',
 'Schizophrenia' : 'teal',
 'TouretteSyndrome' : 'teal',
 'UnipolarDepression' : 'teal'}
row_colors = df.columns.unique().map(lut)

#Genes v diseases with risk scores

df = pd.read_csv('/Users/melis/Documents/Gene-targets/Disease_gene_associations_OpenTargets/risk-score-all-diseases.csv')
df = df.drop_duplicates(subset = ['symbol'])
df = df.set_index('symbol')
x0, _y0, _w, _h = g.cbar_pos

g = sns.clustermap(df, col_colors = row_colors,yticklabels=False)
g.ax_cbar.set_position([x0, 0.02, g.ax_row_dendrogram.get_position().width/6, 0.1])
g.ax_cbar.set_title('genetic risk')
plt.savefig('plots/01-clustering_genes_and_diseases/clustermap_genetic_association.pdf')
plt.savefig('plots/01-clustering_genes_and_diseases/clustermap_genetic_association.png', dpi = 300)
plt.show()

#Total genes overlapping across diseases

df = pd.read_csv('/Users/melis/Documents/Gene-targets/Disease_gene_associations_OpenTargets/hot-encoded-diseases.csv')
df = df.drop_duplicates(subset = ['genes'])
df = df.set_index('genes')

sns.clustermap(df.T.dot(df), vmax = 400, row_colors=row_colors, col_colors = row_colors)
plt.savefig('plots/01-clustering_genes_and_diseases/clustermap_disease_overlaps.png', dpi = 300)
plt.show()

df.T.dot(df).to_csv('processed_data/01-clustering_genes_and_diseases/counts_overlapping_genes.csv')