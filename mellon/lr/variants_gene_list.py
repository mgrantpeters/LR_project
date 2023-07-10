def total_variants_gene_list(genelist):
    '''
    This will calculate the median number of transcript variants for a given list of genes.

    Input:
    genelist - list of genes retrieving variant number for

    Output:
    median number of transcript variants for a gene list
    '''
    import pyensembl
    ensembl = pyensembl.EnsemblRelease(release=109)
    total_transcripts = 0
    for eachgene in genelist:
        total_transcripts += len(ensembl.transcript_ids_of_gene_name(eachgene))
    return total_transcripts/len(genelist)
