#!/usr/bin/env python

'''
This script assembles the gene-risk information for multiple disorders into a one-hot encoded csv. 
INPUTS: 
    - Data downloaded from OpenTargets per disease (.tsv)
PROCESS:
    - Filters the gene lists based on the given risk threshold
    - Checks following this that there are at least 10 genes per disease (otherwise disease is excluded). 
    - Converts into one-hot encoded information (1 gene is associated with disease, 0 gene is not associated with disease)
    - Saves hot one encoded csv 
    - Repeat this process for all given risk thresholds
OUTPUTS: 
    - One-hot-encoded csv with genes as rows and diseases as columns across multiple risk thresholds

'''

import numpy as np
import pandas as pd
import os, glob
from sklearn.preprocessing import OneHotEncoder

path = glob.glob("raw_data/*.tsv")
thr = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
for m in thr:
    for n in range (0,len(path)):
        dis = pd.read_csv(path[n], delimiter = '\t')
        dis = dis[['symbol', 'objectObject']]  # Take gene symbol and column containing genetic association score
        dis = dis[dis['objectObject']>m]  # Filter based on chosen risk threshold 
        disease = glob.glob("raw_data/*.tsv")[n].split("\\")[1].split("_")[0]  # get disease name
        print(disease)
        if np.shape(dis)[0]>10:
                dis['score']=int(1)
                dis = dis.drop(columns=['objectObject']).set_index(['symbol'])

                #creating instance of one-hot-encoder
                encoder = OneHotEncoder(handle_unknown='ignore')
                encoder_df = pd.DataFrame(encoder.fit_transform(dis).toarray())
                encoder_df[disease] = encoder_df[0]
                encoder_df = encoder_df.drop(columns=[0])
                encoder_df['genes'] = list(dis.index)
                encoder_df.set_index('genes', inplace=True)
                if n==0:
                    alldis = encoder_df
                else:
                    alldis = pd.merge(alldis, encoder_df, how='outer', left_index=True, right_index=True).fillna(0)
        else:
            print(disease, " doesn't have enough genes with high genetic scores.")

    alldis.to_csv('processed_data/hot-encoded-diseases_%f.csv' % m)