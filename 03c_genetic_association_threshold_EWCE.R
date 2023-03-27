library("EWCE")
library("biomaRt")
library("ggplot2")
library("grid")
library("dplyr")
set.seed(1234)
setwd("/Users/melis/Documents/SoniaGR_RBP_EWCE/code/")

load("C:/Users/melis/Downloads/ctd_aibsMultipleCrtxSmrtSeq.rda")  #Load CTD data
n = 10000

data = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.1.csv")

#################------------------- Telomere Length
#Load target gene list

list = data[data$Splicing.regulation>0,1]

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

full_results$results$q = p.adjust(full_results$results$p, method = "fdr", n = length(full_results$results$p))
