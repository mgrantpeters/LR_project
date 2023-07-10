#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy

###--------------------------------------------LOAD DATA

#--Load tangram deconvoluted data and combine IF and nonIF data
adata = sc.read("../LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/nonIF_c2l_anndata_combined.h5ad")
adata.obs['key'] = adata.obs.index+'_'+adata.obs['sample_id'].astype(str)
adata.obs['bayesSpace_harmony_9'] = adata.obs['bayesSpace_harmony_9'].astype('category')
adata.obs['bayesSpace_harmony_9'] = adata.obs['bayesSpace_harmony_9'].cat.rename_categories({1: 'SpD1 ~ L1', 
                                                     2: 'SpD2 ~ L1', 
                                                     3 : 'SpD3 ~ L2', 
                                                     4 : 'SpD4 ~ L5', 
                                                     5 : 'SpD5 ~ L3', 
                                                     6 : 'SpD6 ~ WM', 
                                                     7 : 'SpD7 ~ L6', 
                                                     8 : 'SpD8 ~ L4', 
                                                     9 : 'SpD9 ~ WM'})

###--------------------------------------------RUN ANALYSIS

reslist = ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7']

for res in range(0,len(reslist)):
    print(reslist[res])
    top = pd.read_csv('processed_data/03-LR_network_visualisation/louvain_largest_cluster_%s.csv' % reslist[res], sep = ',', index_col = 0)
    for n in range(0,1000):
        sc.tl.score_genes(adata, list(adata.var[adata.var['gene_name'].isin(list(top['0']))].index), ctrl_size = len(list(top['0'])))
        df = adata.obs['score']
        if (n == 0):
            all_scores = df
        else: 
            all_scores = pd.concat([all_scores, df], axis = 1)
    adata.obs['score'] = all_scores.median(axis = 1)
    if (res==0):
        temp = adata.obs[['bayesSpace_harmony_2', 'bayesSpace_harmony_9','bayesSpace_harmony_16','bayesSpace_harmony_28', 'score']]
        temp['resolution'] = reslist[res]
        sd = temp
    else:
        temp = adata.obs[['bayesSpace_harmony_2', 'bayesSpace_harmony_9','bayesSpace_harmony_16','bayesSpace_harmony_28', 'score']]
        temp['resolution'] = reslist[res]
        sd = pd.concat([sd,temp], axis = 0)
    sns.set_theme(style="whitegrid", palette= sns.color_palette("husl", 9))
    plt.figure(figsize=(10,7), dpi = 300)
    sns.boxplot(data=sd, x="bayesSpace_harmony_9", y="score", showfliers=False)
    plt.title('enrichment score per spatial domain - resolution %s' % reslist[res])
    plt.xticks(rotation=90)
    plt.savefig('plots/03-LR_network_visualisation/03e_network_spatial_analysis/res_%s_bayespace_harmony_9_across_resolutions.pdf' % reslist[res], bbox_inches = 'tight')
    plt.show()

    pvals = []
    for cluster in sd['bayesSpace_harmony_9'].unique().sort_values():
        subset_res = sd[sd['resolution']==reslist[res]]
        subset = subset_res[subset_res['bayesSpace_harmony_9']==cluster]
        pvals.append((subset['score']<=0).sum()/len(subset['score']))

    print(pd.DataFrame(pvals, index = sd['bayesSpace_harmony_9'].unique().sort_values()))

    pd.DataFrame(pvals, index = sd['bayesSpace_harmony_9'].unique().sort_values()).to_csv('processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/res_%s_pvalues_scores_per_cluster.csv' % reslist[res])
    sd.to_csv('processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/scores_per_spatial_domain_all_resolutions.csv')
    
sns.set_theme(style="whitegrid", palette= sns.color_palette("husl", 9))
plt.figure(figsize=(10,7), dpi = 300)
sns.lineplot(data=sd.reset_index(), x="resolution", y="score", hue = "bayesSpace_harmony_9")
plt.title('Network score for spatial domains across genetic risk thresholds')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('plots/03-LR_network_visualisation/03e_network_spatial_analysis/bayespace_harmony_9_across_resolutions.pdf', bbox_inches = 'tight')
plt.show()
