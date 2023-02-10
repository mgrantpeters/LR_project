def find_hits(target, db):                                                                      
    '''
    
    Function for identifying genes that are targets or receptors (tgts) within a given database (db).
    
    Input: 
    target - Pandas dataframe where the targets of interest are in the index; name of index must be 'gene'
    db - Reference ligand receptor database; expects columns to match omnipath database

    Output: 
    final_s - List of ligands in the given target list (termed 'source' by omnipath)
    perc_Src - percentage of genes that act as source
    final_t - List of receptors in the given target list (termed 'targets' by omnipath)
    percTgt - percentage of genes that act as targets
    Sum - Variable with information regarding whether a gene acts as ligand or receptor (1) or both (2)
    
    '''
    import pandas as pd 

    hits = db['genesymbol_intercell_source'].isin(list(target.index))                           #Out of targets, which act as sources (aka: ligands)? 
    hits = hits.to_frame()
    hits['gene'] = list(db['genesymbol_intercell_source'])
    hits = hits[hits['genesymbol_intercell_source']==True]
    
    print("Source hits: ", len(hits['gene'].unique()))
    if (len(target.index.unique())>0):
        print("% of source hits in total genes of interest: ", len(hits['gene'].unique())/len(target.index.unique())*100)
        percSrc = len(hits['gene'].unique())/len(target.index.unique())
    else:
        percSrc = 0

    final_s = pd.merge(hits, db[['transmitter_intercell_source', 'receiver_intercell_source',
        'secreted_intercell_source',
        'plasma_membrane_transmembrane_intercell_source',
        'plasma_membrane_peripheral_intercell_source']], how="inner", left_index = True, right_index = True)
    final_s['source']=1
    final_s[['gene', 'source']]

    hits = db['genesymbol_intercell_target'].isin(list(target.index))
    hits = hits.to_frame()
    hits['gene'] = list(db['genesymbol_intercell_target'])
    hits = hits[hits['genesymbol_intercell_target']==True]
    
    print("Target hits: ", len(hits['gene'].unique()))
    if (len(target.index.unique())>0):
        print("% of target hits in total genes of interest: ", len(hits['gene'].unique())/len(target.index.unique())*100)
        percTgt = len(hits['gene'].unique())/len(target.index.unique())*100
    else:
        percTgt=0

    final_t = pd.merge(hits, db[['transmitter_intercell_target', 'receiver_intercell_target',
        'secreted_intercell_target',
        'plasma_membrane_transmembrane_intercell_target',
        'plasma_membrane_peripheral_intercell_target']], how="inner", left_index = True, right_index = True)
    final_t['target']=2
    final_t[['gene', 'target']]
    
    target['Empty'] = 0
    all = target['Empty'].reset_index()
    all.rename(columns = {'index':'gene'}, inplace = True)
    summary_LR = pd.merge(final_s[['gene', 'source','secreted_intercell_source',
       'plasma_membrane_transmembrane_intercell_source',
       'plasma_membrane_peripheral_intercell_source']], final_t[['gene', 'target', 'secreted_intercell_target',
       'plasma_membrane_transmembrane_intercell_target',
       'plasma_membrane_peripheral_intercell_target']], how="outer", on=['gene','gene']).fillna(0)
    summary_LR = pd.merge(summary_LR.drop_duplicates(subset = ["gene"]), all, how="outer", on = ['gene','gene']).fillna(0)
    summary_LR['source+target']=summary_LR['source']+summary_LR['target']
    summary_LR['source+target'].astype("category")
    summary_LR['source+target'].replace({
    0.0: 'None', 
    1.0: 'Source', 
    2.0: 'Target', 
    3.0: 'Source+Target'}, inplace = True)
    
    return(final_s, percSrc, final_t, percTgt, summary_LR.drop_duplicates(subset = ["gene"]))