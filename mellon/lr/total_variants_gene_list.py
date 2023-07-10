def total_variants_gene_list(genelist):
    '''
    Retrieves from ENSEMBL release 109 the number of known transcripts for each gene in a given list

    Input:
    genelist - list of genes you want the median transcripts for

    Output:
    median of total transcripts in the given list
    '''

    import random
    import mellon as ml
    import pyensembl
    import numpy as np

    ensembl = pyensembl.EnsemblRelease(release=109)
    total_transcripts = []
    for eachgene in genelist:
        try:
            total_transcripts.append(len(ensembl.transcript_ids_of_gene_name(eachgene)))
        except:
            print(eachgene, ' is not available in ENSEMBL. Excluding this gene from analysis.')
            continue
    return np.median(total_transcripts)
