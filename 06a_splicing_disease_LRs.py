#!/usr/bin/env python

'''
This script fetches from ENSEMBL (release 109) the number of transcript variants of the LRs of interest per disease. 
It then compares the number of transcript variants in that disease to bootstrapped data of LRs expressed in the brain of the same size,
calculating a p-value for this. 

INPUTS: 
    - GTEx median gene expression (v 7)
    - OmnipathDB
    - ENSEMBL 109
    - OpenTargets risk genes .tsv
PROCESS:
    - 'universe' list: Fetch all genes expressed in the brain (GTEx) and filter to only include those which act as LRs
    - calculate median number of transcript variants for LRs per disease per risk threshold
    - Use universe to bootstrap lists of LRs expressed in brain, then calculate the median transcript variants for each list
    - calculate pval for inquiry data in relation to bootstrapped data

OUTPUTS: 
    - .csv file with information about disease, median transcript variants, pvalue and risk threshold

'''


import numpy as np
import pandas as pd
import mygene
import omnipath as op
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import mellon as ml

gtex_link = 'https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz'
exp = pd.read_csv(gtex_link, sep='\t', index_col='gene_id', skiprows=2)
exp_cns = exp.loc[:, ['Brain - Anterior cingulate cortex (BA24)',
       'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
       'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
       'Brain - Hippocampus', 'Brain - Hypothalamus',
       'Brain - Nucleus accumbens (basal ganglia)',
       'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
       'Brain - Substantia nigra']]
exp_cns = exp_cns.loc[(exp_cns > 0).any(axis=1)]

# import omnipath db
db = op.interactions.import_intercell_network(transmitter_params = {"categories":"ligand"}, receiver_params = {"categories": "receptor"})
db = db[np.logical_not(db['genesymbol_intercell_source'].str.startswith('HLA'))]
db = db[np.logical_not(db['genesymbol_intercell_target'].str.startswith('HLA'))]
db = db[~db['genesymbol_intercell_target'].astype(str).str.startswith('COMPLEX')]
db = db[~db['genesymbol_intercell_source'].astype(str).str.startswith('COMPLEX')]

# convert gtex gene names from ENSEMBL to gene symbols
mg = mygene.MyGeneInfo()
ensembl_gtex = list(np.unique(pd.DataFrame(list(exp_cns.index.str.split('.')))[0]))
symbols_gtex = mg.querymany(ensembl_gtex, scopes='ensembl.gene', fields='symbol', species='human')
symbols_gtex = pd.DataFrame(symbols_gtex)['symbol']

# Keep only ones that act as LRs
symbols_gtex_format = symbols_gtex.rename('gene').to_frame().set_index('gene')
all_source, pctsrc, all_targets, pcttgt, summary_LR = ml.lr.find_hits(symbols_gtex_format,db) 
universe = list(summary_LR[summary_LR['source+target']!='None']['gene'])

rng = [0.4, 0.7, 0.1]

for m in rng:
    #fetch disease risk list path
    path = glob.glob("/Users/melis/Documents/Gene-targets/Disease_gene_associations_OpenTargets/*.tsv")
    dis_list = []
    #iterate through each disease
    for n in range (0,len(path)):
        disease = glob.glob("/Users/melis/Documents/Gene-targets/Disease_gene_associations_OpenTargets/*.tsv")[n].split("\\")[1].split("_")[0]
        dis = pd.read_csv(path[n], delimiter = '\t', index_col=0)
        dis = dis[(dis['objectObject']>m)]
        # keep only diseases that, at this threshold, have 10 or more genetic associations
        if np.shape(dis)[0]>10:
            dis = pd.read_csv(path[n], delimiter = '\t', index_col=0)
            dis = dis[(dis['objectObject']>m)]
            if np.shape(dis)[0]>10:
                df = dis[(dis['objectObject']>m)]
                df['gene']=df.index
                df.reset_index(inplace = True)
                df = df.drop(columns=['symbol']).set_index(['gene'])
                # calculate the occurance of Ls and Rs in our disease data
                all_source, pctsrc, all_targets, pcttgt, summary_LR = ml.lr.find_hits(df,db)
                # keep list of unique genes that act as L and/or R, calculate number of transcript variants 
                LRs = list(summary_LR[summary_LR['source+target']!='None']['gene'])
                disease_med_var = ml.lr.total_variants_gene_list(LRs)
                # bootstrap random list of LRs expressed in the brain (minus targets), calculate number of transcript variants and store these in a list 
                vars = ml.lr.bootstrap_genes_return_variants(len(LRs), 10000, len(set(universe) - set(LRs)))

                plt.hist(vars, bins=20)
                plt.show()
                
                #We calculate the pval
                pval = np.sum(pd.DataFrame(vars)>disease_med_var)/len(vars)

                #vars_all_diseases.append(disease_med_var)
            dis_list.append(disease)

            if (n==0):
                all_disease_var = np.asarray(disease_med_var).reshape(1,-1)
                all_pval = pval
                
            else:    
                all_disease_var = np.concatenate((all_disease_var, np.asarray(disease_med_var).reshape(1,-1)), axis = 0)
                all_pval = np.concatenate((all_pval, pval), axis = 0)
    threshold = np.zeros_like(all_pval)
    threshold[:] = m

    # Store in a dataframe the disease name, number of transcripts (total or median per gene?) and the p-value
    
    results_df = pd.DataFrame(np.concatenate((np.asarray(dis_list).reshape(1,-1), all_disease_var.reshape(1,-1), all_pval.reshape(1,-1), threshold.reshape(1,-1)), axis = 0).T, columns = ['Disease', 'Median transcript variants', 'pval', 'Threshold'])
    if m == rng[0]:
        all_results_df = results_df

    else:
        all_results_df = pd.concat([all_results_df, results_df], axis = 0)

all_results_df.to_csv('processed_data/06a_splicing_disease_LRs/significance_splicing_disease.csv', index = False)
