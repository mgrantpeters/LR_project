#!/usr/bin/env bash

'''
This script is designed to run in an automated manner the entire analysis of this repository. 
To find out more about the python and R requirements and data requirements, read our README file. 

'''

# Initial hypothesis testing: are LRs likely to be disrupted in neurological and neuropsychiatric disease?

python 00-hot_one_encode_disease_genes.py       # Hot encode
python 01-clustering_genes_and_diseases.py      # Cluster genes v diseases, counts of directly overlaping genes
python 01a_LR_occurrance_across_diseases.py     # Do LRs occur more frequently in LRs associated with disease than other brain-expressed LRs?
Rscript 01b_plot_LR_per_disease.R
python 06a_splicing_disease_LRs.py              # Are disease LR interactions more vulnerable to mis-splicing? 
Rscript 06b_plot_splicing_results.R
python 02-clustering_with_LR_partners.py        # If we cluster LR pairs v diseases, are diseases more similar to each other?
python 03_LR_network_distance_diseases.py       # Calculating the distance between LR networks of each disease

#####---------- Still to make .py scripts for:
# The cross-disease network
# Characterising risk score thresholds
python 03a_genetic_association_test.py          # Characterise thresholds, generate a network per threshold
python 03d_network_stats_per_threshold.py       # Are the networks bigger or smaller than we expect? Do they have better or worse connectivity?
python 03f_characterising_diseases_in_networks.py # LRs associated to each disease per threshold

# Cell type enrichment
Rscript 03c_a_genetic_association_threshold_EWCE_annotation1.R 

# GO pathway enrichment
Rscript 03h_GO_enrich_across_thresholds.R

# Spatial characterisation
python 03e_network_spatial_analysis.py          # Which spatial domains are enriched for network?
python 03-1e_visualising_spatial_boxplots.py    # Visualise boxplots of script 03e
python 03h_spatial_colocalisation.py
python 03i_histograms_detectability_network_genes_spatial.py
python 03j_spatial_plots_scores.py

# Visualisations
python 03g_comparing_disease_networks.py        # Visualise two diseases in the same network
python 03i_hotspot_visalisation.py              # May remove this