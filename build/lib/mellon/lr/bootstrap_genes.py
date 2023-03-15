def bootstrap_genes(n, iter, genelist, db):
    '''
    Function for creating a bootstrapped gene list with n genes drawn randomly from genelist.
    This will then identify how many are ligands or receptors from a database db.

    Input:
    n - number of genes
    iter - number of iterations
    genelist - gene pool to randomly bootstrap genes from
    db - database to detect Ls and Rs from

    Output:
    src - bootstrapped occurance of sources
    tgt - bootstrapped occurance of ligands
    '''

    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    import mellon as ml
    src = []
    tgt = []

    for i in range (0,iter):
        rand = random.sample(range(0,len(genelist)), n)                                         # Generate n random numbers for gene selection
        source, percSrc, targets, percTgt, sum = ml.lr.find_hits(genelist.iloc[rand], db)      # Select genes
        src.append(percSrc)
        tgt.append(percTgt)

    sns.histplot(data = src, bins = 20)
    plt.show()

    sns.histplot(data = tgt, bins = 20)
    plt.show()

    return src, tgt
