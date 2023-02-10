def topocurring(targets, hits):                                                                 #Function for identifying top occuring genes in a set of groups (e.g.: diseases, databases)
    '''
    Function for identifying genes that are targets or receptors (tgts) within a given database (db).
    
    Input: 
    targets - pandas dataframe with index of genes and a column per disease of interest
    hits - source or target output from find_hits to filter for sources/targets

    Output: 
    top10S - Returns a pandas dataframe with list of genes, number of diseases it is a risk for and the list of diseases it is risk for
    '''
    import numpy as np

    top10S = targets.loc[list(hits['gene'].unique())].sum(axis = 1).sort_values(ascending = False).to_frame(name = 'Occurances')
    top10S['Diseases'] = 'A'
    for n in range (0, len(list(top10S.index))):
        info = targets.loc[top10S.index[n]]
        info = info[info>0].to_frame()
        l = []
        for i in range (0,len(info.index)):
            l.append(info.index[i].split("_")[0])
        top10S['Diseases'][n] = ", ".join(np.unique(l))
    return top10S