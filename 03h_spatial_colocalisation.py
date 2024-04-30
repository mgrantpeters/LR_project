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
import pickle

###--------------------------------------------LOAD DATA
adata1 = sc.read("../LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/nonIF_c2l_anndata_combined_30samples.h5ad")
adata2 = sc.read("../LIBD_LR/grantpm/DLPFC_Visium_LIBD/processed-data/MGP_analysis/IF_c2l_anndata_combined.h5ad")

adata = adata1.concatenate(adata2, join = 'inner')
if (adata.var['gene_name-1'].astype(str).equals(adata.var['gene_name-0'].astype(str))==True):
    adata.var['gene_name'] = adata.var['gene_name-1']
    adata.var.drop(columns = ['gene_name-1', 'gene_name-0'], inplace = True)

def coloc_LR_list(original, gene_list, res):
    adata = copy.copy(original)
    counter = 0
    for n in range(0,len(gene_list)):
        gene = adata.var['gene_name'].eq(gene_list[n])
        if (gene.sum()>0):
            counter+=1
            if (counter==1):
                df = (adata.to_df()[adata.var[gene].index]>0).astype(int)
            else: 
                df[adata.var[gene].index[0]] = (adata.to_df()[adata.var[gene].index]>0).astype(int)
                
    df['sum'] = df.sum(axis=1)
    print(counter)
    print(counter/len(gene_list)*100)
    df['sum'].hist(bins = len(gene_list))
    df.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_dataframe_counts_expression_network.csv" % (res))
    return df

def hot1_encode_cell_counts(adata):
    probs = adata.obs[['Astro', 'EndoMural', 'Excit_L2_3', 'Excit_L3',
       'Excit_L3_4_5', 'Excit_L4', 'Excit_L5', 'Excit_L5_6', 'Excit_L6',
       'Inhib', 'Micro', 'OPC', 'Oligo']]

    cells = []

    for n in range (0,len(probs.columns)):
        cells.append(probs.columns[n].replace('prop_',''))
    probs.columns = cells

    print("Building adjacency matrix for top 3 cells per spot...")
    hot1 = copy.copy(probs)
    hot1.iloc[:,:] = 0
    for m in range (0,np.shape(probs)[0]):
        hot1.iloc[m,:][probs.T[probs.index[m]].nlargest(3).index]=1
    return hot1


def network_top3(adata, tgt, pvals):
    probs = adata.obs[['Astro', 'EndoMural', 'Excit_L2_3', 'Excit_L3',
       'Excit_L3_4_5', 'Excit_L4', 'Excit_L5', 'Excit_L5_6', 'Excit_L6',
       'Inhib', 'Micro', 'OPC', 'Oligo']]

    cells = []

    for n in range (0,len(probs.columns)):
        cells.append(probs.columns[n].replace('prop_',''))
    probs.columns = cells

    cum_prob = probs.max(axis = 1).to_frame()
    cum_prob.rename(columns={0: "Probability"}, inplace = True)
    cum_prob['n'] = 1
    cum_prob = cum_prob.to_numpy()
    print("Calculating cumulative probability...")
    for n in range (1,20):
        for m in range (0,np.shape(probs)[0]):
            add = [probs.T[probs.index[m]].nlargest(n).sum(), n]
            add = np.asarray(add).reshape(1,-1)
            cum_prob = np.concatenate((cum_prob, add), axis = 0)

    print("Done")        
    df = pd.DataFrame(cum_prob, 
                 columns=['Probability', 
                          'n'])

    ax = sns.violinplot(x="n", y="Probability",
                        data=df, palette="Set2")
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_cumulative_probability.pdf" % (tgt), dpi = 150, bbox_inches = 'tight')
    plt.show()
    plt.clf()

    #for n in range (1,20):
    #    print(df[df['n']==n].mean())

    print("Building adjacency matrix for top 3 cells per spot...")
    hot1 = copy.copy(probs)
    hot1.iloc[:,:] = 0
    for m in range (0,np.shape(probs)[0]):
        hot1.iloc[m,:][probs.T[probs.index[m]].nlargest(3).index]=1

    adj = hot1.T.dot(hot1)
    np.fill_diagonal(adj.values, 0)
    adjmat_notnorm = adj
    adjmat = (adj/adj.sum().sum())*100
    print("Done.")
    #print(adjmat)

    G = nx.from_numpy_matrix(adjmat.values)  
    edges = dict(zip(list(range(0,len(adjmat.columns))), adjmat.columns))
    
    nodes_list = []
    listing = [list(i) for i in list(G.edges)]
    
    for m in range(0,len(list(G.nodes))):
        counter = 0
        for n in range(0,len(listing)):
            if m in listing[n]:
                counter +=1
            #print(counter)
        #if counter == 0:
            #del color_map[m]           
        if counter!= 0:
            nodes_list.append(m)     
            
    #del(color_map)
    color_map = []
    node_size = []
    for i in range(0,len(nodes_list)):
        if nodes_list[i] == 0:
            color_map.append('yellowgreen')
        if nodes_list[i] == 1:
            color_map.append('darkorange')
        if nodes_list[i] == 2 or nodes_list[i] == 3 or nodes_list[i] == 4 or nodes_list[i] == 5 or nodes_list[i] == 6 or nodes_list[i] == 7 or nodes_list[i] == 8:
            color_map.append('plum')
        if nodes_list[i] == 9:
            color_map.append('thistle')
        #if (nodes_list[i] == 4):
        #    color_map.append('darkgreen')
        if nodes_list[i] == 10:
            color_map.append('indianred')
        if nodes_list[i] == 11:
            color_map.append('skyblue')
        if nodes_list[i] == 12:
            color_map.append('cornflowerblue')
    

    for i in range(0,len(cells)):
        if ((adjmat.loc[[cells[i]], :].sum().sum()+adjmat.loc[:, [cells[i]]].sum())>0)[0]:
            node_size.append(250*(adjmat.loc[[cells[i]], :].sum().sum()+adjmat.loc[:, [cells[i]]].sum()))
    print(node_size)
    print("Plotting...")
    esmallest = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] < 0]
    esmall = [(u, v) for (u, v, d) in G.edges(data=True) if 0 <= d["weight"] < 1]
    emedium = [(u, v) for (u, v, d) in G.edges(data=True) if 1 <= d["weight"] < 2]
    elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] >= 2]
    
    if (len(pvals.columns)!=0):
        pvalT = pvals.T
        print(pvals.T)
        pvalT.loc[pvalT[0] <= 0.05, 'alpha'] = 1
        pvalT.loc[pvalT[0] > 0.05, 'alpha'] = 0.3

    plt.figure(figsize=(20, 18))
    pos = nx.spring_layout(G, seed=7)  # positions for all nodes - seed for reproducibility
    if (len(pvals.columns)!=0):
        nx.draw_networkx_nodes(G, pos, nodelist = nodes_list, node_color = color_map, node_size=node_size, alpha = list(pvalT['alpha']))
    else:
        nx.draw_networkx_nodes(G, pos, nodelist = nodes_list, node_color = color_map, node_size=node_size)
    nx.draw_networkx_edges(G, pos, edgelist=esmallest, width=0.5, alpha=0.05, edge_color="black", style="dashed")
    nx.draw_networkx_edges(G, pos, edgelist=esmall, width=2, alpha=0.5, edge_color="black", style="dashed")
    nx.draw_networkx_edges(G, pos, edgelist=emedium, width=6, alpha=0.7, edge_color="black")
    nx.draw_networkx_edges(G, pos, edgelist=elarge, width=9, edge_color="black")
    nx.draw_networkx_labels(G, pos, font_size=18, labels = edges, verticalalignment = 'top', font_family="sans-serif")
    #nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")
    #nx.draw(G)
    ax = plt.gca()
    ax.margins(0.1)
    #plt.axis("off")
    plt.tight_layout()
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_network.pdf" % (tgt), dpi = 150, bbox_inches = 'tight')
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_network.png" % (tgt), dpi = 300)
    plt.show()
    plt.clf()
    
    pickle.dump(G, open('processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_network.pickle' % (tgt), 'wb'))
    pd.DataFrame(color_map).to_csv('processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_network_colormap.csv' % (tgt))
    pd.DataFrame(list(node_size)).to_csv('processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_network_nodesize.csv' % (tgt))


    #sns.clustermap(full_adjmat, cmap='Spectral_r', row_colors=color_map, col_colors=color_map, row_cluster=False, col_cluster = False)
    #plt.savefig("Network_analysis/%s_localisation_heatmap.png" % (groups[g]), dpi = 150)
    #plt.show()
    
    print("Done.")

    print("PageRank is:")
    G_notnorm = nx.from_numpy_matrix(adjmat_notnorm.values)  
    pgrank = pd.DataFrame.from_dict(nx.pagerank(G_notnorm, alpha=0.85, weight = 'weight'), orient='index')
    pgrank.index = cells
    pgrank.sort_values(by=[0], ascending = False, inplace = True)
    print(pgrank)
    pgrank.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_PageRank_alpha=0.85.csv" % (tgt))
    
    sns.heatmap(pd.DataFrame(np.tril(adjmat), index = adjmat.index, columns = adjmat.columns), cmap = 'inferno_r')
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_heatmap.pdf" % (tgt), dpi = 150, bbox_inches = 'tight')
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_heatmap.png" % (tgt), dpi = 300)
    plt.show()
    plt.clf()
    adjmat_notnorm.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_adjacencymatrix_notnormalised.csv" % (tgt))
    adjmat.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_c2l_3cells_adjacencymatrix.csv" % (tgt))
    return adjmat, hot1.sum()

###--------------------------------------------RUN ANALYSIS FOR TOP 3 CELLS

res = ['0.7', '0.55', '0.4']
iter = 1000
qtl = 0.98
for n in range(0,np.shape(res)[0]):
    top = pd.read_csv('processed_data/03-LR_network_visualisation/louvain_largest_cluster_%s.csv' % res[n], sep = ',', index_col = 0)
    pvals=[]
    df = coloc_LR_list(adata, list(top['0']), res[n])

    interest = adata[adata.obs.index.isin(df['sum'][df['sum']>df['sum'].quantile(q=qtl)].index)]
    interest.obs['bayesSpace_harmony_9'] = interest.obs['bayesSpace_harmony_9'].cat.rename_categories({1: 'SpD1 ~L1', 
                                                         2 : 'SpD2 ~L1', 
                                                         3 : 'SpD3 ~ L2', 
                                                         4 : 'SpD4 ~ L5', 
                                                         5 : 'SpD5 ~ L3', 
                                                         6 : 'SpD6 ~ WM', 
                                                         7 : 'SpD7 ~ L6', 
                                                         8 : 'SpD8 ~ L4', 
                                                         9 : 'SpD9 ~ WM'})
    data = pd.DataFrame(interest.obs[['sample_id', 'bayesSpace_harmony_9']].value_counts(normalize=True)).reset_index().rename(columns={0 : '% spots'})
    data = data.dropna()
    sns.set_theme(style="white", palette= sns.color_palette("husl", 9))
    sns.boxplot(data=data, x="bayesSpace_harmony_9", y="% spots", showfliers=False, order = ['SpD1 ~ L1', 'SpD2 ~ L1', 'SpD3 ~ L2','SpD5 ~ L3','SpD8 ~ L4', 'SpD4 ~ L5', 'SpD7 ~ L6', 'SpD6 ~ WM', 'SpD9 ~ WM'])
    plt.xticks(rotation=90)
    plt.savefig("plots/03-LR_network_visualisation/03h_spatial_colocalisation/%s_boxplot_c2l_top2pct_perlayer.png" % (res[n]), bbox_inches = 'tight')
    plt.show()
    data.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_data_c2l_top2pct_perlayer.csv" % (res[n]))
    LRadjmat, cell_counts = network_top3(interest, 'LR_'+res[n], pd.DataFrame([]))
    for i in range(0,iter):
        number_spots = len(interest.obs.index)
        rand_spots = random.sample(list(adata.obs.index), number_spots)
        hot1 = hot1_encode_cell_counts(adata[adata.obs.index.isin(rand_spots)])
        if (i == 0):
            bootstrapped = hot1.sum()
        else:
            bootstrapped = pd.concat([bootstrapped, hot1.sum()], axis = 1)
    bootstrapped = bootstrapped.T
    counts = pd.DataFrame(cell_counts).T
    for celltype in hot1.columns:
        pvals.append((1+((bootstrapped[celltype]>=counts[celltype][0]).sum()))/(iter+1))
    final_pvals = pd.DataFrame(np.asarray(pvals).reshape(1,-1), columns = hot1.columns)
    final_pvals.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_bootstrapped_cell_type_pvals.csv" % (res[n]))
    bootstrapped.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_bootstrappeddata_cell_type.csv" % (res[n]))
    counts.to_csv("processed_data/03-LR_network_visualisation/03h_spatial_colocalisation/%s_realcounts_cell_type.csv" % (res[n]))
    
    LRadjmat, cell_counts = network_top3(interest, 'LR_'+res[n], final_pvals)