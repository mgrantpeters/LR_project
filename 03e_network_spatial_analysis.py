#!/usr/bin/env python

###--------------------------------------------LOAD LIBRARIES
import scanpy as sc
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
import random
from scipy.stats import mannwhitneyu

###--------------------------------------------LOAD DATA

#--Load tangram deconvoluted data and combine IF and nonIF data
adata = sc.read("../LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/nonIF_c2l_anndata_combined_30samples.h5ad")
adata.obs['key'] = adata.obs.index+'_'+adata.obs['sample_id'].astype(str)
adata.obs['bayesSpace_harmony_9'] = adata.obs['bayesSpace_harmony_9'].astype('category')
adata.obs['bayesSpace_harmony_9'] = adata.obs['bayesSpace_harmony_9'].cat.rename_categories({1: 'L1 ~ SpD1', 
                                                     2 : 'L1 ~ SpD2', 
                                                     3 : 'L2 ~ SpD3', 
                                                     4 : 'L5 ~ SpD4', 
                                                     5 : 'L3 ~ SpD5', 
                                                     6 : 'WM ~ SpD6', 
                                                     7 : 'L6 ~ SpD7', 
                                                     8 : 'L4 ~ SpD8', 
                                                     9 : 'WM ~ SpD9'})
adata.obs["bayesSpace_harmony_9"] = (adata.obs["bayesSpace_harmony_9"].cat.reorder_categories(sorted(adata.obs["bayesSpace_harmony_9"].cat.categories)))

###--------------------------------------------RUN ANALYSIS

reslist = ['0.4', '0.55', '0.7']

for res in range(0,len(reslist)):
    print(reslist[res])
    top = pd.read_csv('processed_data/03-LR_network_visualisation/louvain_largest_cluster_%s.csv' % reslist[res], sep = ',', index_col = 0)
    for n in range(0,1000):
        #Changed the scores to ensure the universe does not resample the tested genes, as was the scanpy default
        genes_of_interest = list(adata.var[adata.var['gene_name'].isin(list(top['0']))].index)
        universe = list(adata.var.index)
        #gene_pool = list(set(universe) - set(set(universe) & set(genes_of_interest)))   #establish gene pool as universe minus target genes
        gene_pool = list(set(universe))
        interest = adata[:, genes_of_interest].to_df().sum(axis = 1)                    #Take the sum of gene expression for genes of interest                     
        r = random.sample(gene_pool, len(genes_of_interest))
        background = adata[:, r].to_df().sum(axis = 1)                                  #Take sum of expression for background genes
        adata.obs['score'] = (interest - background)/adata.obs['sum_umi']               #Normalise expression by total UMIs captured per spot
        adata.obs['score'] = (interest - background)/adata.obs['sum_umi']
        adata.obs['genes_of_interest'] = interest
        adata.obs['background'] = background
        df1 = adata.obs['score']
        df2 = adata.obs['genes_of_interest']
        df3 = adata.obs['background']
        if (n == 0):
            all_scores = df1
            all_geneofinterest = df2
            all_background = df3
        else: 
            all_scores = pd.concat([all_scores, df1], axis = 1)
            all_genesofinterest = pd.concat([all_geneofinterest, df2], axis = 1)
            all_background = pd.concat([all_background, df3], axis = 1)
            
    adata.obs['score_final'] = all_scores.median(axis = 1)                               #take median of all scores. This value was normalised by UMIs as scores across layers are going to be compared                  
    adata.obs['genes_of_interest_final'] = all_genesofinterest.median(axis = 1)          #Store median genes of interest
    adata.obs['background_final'] = all_background.median(axis = 1)                      # Store median background genes

    adata.obs[['score_final','genes_of_interest_final','background_final', 'bayesSpace_harmony_9']].to_csv('plots/03-LR_network_visualisation/03e_network_spatial_analysis/res_%s_scores_genes_interest.csv' % reslist[res])

    if (res==0):
        temp = adata.obs[['bayesSpace_harmony_2', 'bayesSpace_harmony_9','bayesSpace_harmony_16','bayesSpace_harmony_28', 'score_final']]  #Bug fixed on 17/01 - this was taking score, not score_final
        temp['resolution'] = reslist[res]
        sd = temp
    else:
        temp = adata.obs[['bayesSpace_harmony_2', 'bayesSpace_harmony_9','bayesSpace_harmony_16','bayesSpace_harmony_28', 'score_final']]
        temp['resolution'] = reslist[res]
        sd = pd.concat([sd,temp], axis = 0)
    sns.set_theme(style="whitegrid", palette= sns.color_palette("husl", 9))
    plt.figure(figsize=(10,7), dpi = 300)
    sns.boxplot(data=sd, x="bayesSpace_harmony_9", y="score_final", showfliers=False)
    plt.title('enrichment score per spatial domain - resolution %s' % reslist[res])
    plt.xticks(rotation=90)
    plt.savefig('plots/03-LR_network_visualisation/03e_network_spatial_analysis/res_%s_bayespace_harmony_9_across_resolutions.pdf' % reslist[res], bbox_inches = 'tight')
    plt.show()

    for cluster in sd['bayesSpace_harmony_9'].unique().sort_values():
        #Added 18/12/23 to visualise new scoring 
        fig, ax = plt.subplots()
        sns.histplot(data=adata.obs[adata.obs['bayesSpace_harmony_9']==cluster], x="genes_of_interest_final", alpha = 0.6)
        sns.histplot(data=adata.obs[adata.obs['bayesSpace_harmony_9']==cluster], x="background_final", color = 'grey', alpha = 0.6)
        ax.set_xlim([-30, 200])
        plt.savefig('processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/histogram_%s_%s_scores_per_cluster.png' % (reslist[res],cluster), bbox_inches = 'tight', dpi=300)
        plt.show()
        
        adata.obs[adata.obs['bayesSpace_harmony_9']==cluster][['genes_of_interest_final','background_final']].to_csv('processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/scores_for_stats_%s_%s.csv' % (reslist[res],cluster))
        
        subset_res = sd[sd['resolution']==reslist[res]]
        subset = subset_res[subset_res['bayesSpace_harmony_9']==cluster]

    sd.to_csv('processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/scores_per_spatial_domain_all_resolutions.csv')
    
sns.set_theme(style="whitegrid", palette= sns.color_palette("husl", 9))
plt.figure(figsize=(10,7), dpi = 300)
sns.lineplot(data=sd.reset_index(), x="resolution", y="score_final", hue = "bayesSpace_harmony_9")
plt.title('Network score for spatial domains across genetic risk thresholds')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('plots/03-LR_network_visualisation/03e_network_spatial_analysis/bayespace_harmony_9_across_resolutions.pdf', bbox_inches = 'tight')
plt.show()