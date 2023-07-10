def bootstrap_genes_return_variants(n, iter, universe):
    '''
    Function for creating a bootstrapped gene list with n genes drawn randomly from genelist.
    This will then calculate the median number of transcript variants.

    Input:
    n - number of genes
    iter - number of iterations
    universe - gene pool to randomly bootstrap genes from

    Output:
    vars - bootstrapped list of median number of transcript variants per iteration

    '''

    import random
    import mellon as ml
    import pandas as pd
    vars = []

    for i in range (0,iter):
        rand = random.sample(range(0,len(universe)), n)                                         # Generate n random numbers for gene selection
        subset_universe = list(pd.DataFrame(universe).iloc[rand][0])
        med_var = ml.lr.total_variants_gene_list(subset_universe)                           # Select genes
        vars.append(med_var)

    return vars
