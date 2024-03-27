library("EWCE")
library("biomaRt")
library("ggplot2")
library("grid")
library("dplyr")
set.seed(1234)
setwd("/Users/melis/Documents/GitHub/LR_project")

load("C:/Users/melis/Downloads/ctd_aibsMultipleCrtxSmrtSeq.rda")  #Load CTD data
n = 10000


#################------------------- 0.7 threshold
data = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.7.csv", row.names = 1)
list = data$X0
#Load target gene list

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length
thisResult1 = data.frame(full_results$results)
thisResult1$MajorCellType = unlist(lapply(strsplit(as.character(thisResult1$CellType),'_'), `[[`, 1))
thisResult1$testList = '0.7 threshold'
FinalResult1 = thisResult1


plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation1_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/annotation1_thr07.csv", row.names=TRUE, quote=FALSE) 

rm(full_results)
rm(plot_list)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 2,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

thisResult = data.frame(full_results$results)
thisResult$MajorCellType = unlist(lapply(strsplit(as.character(thisResult$CellType),'_'), `[[`, 1))
thisResult$testList = '0.7 threshold'
FinalResult = thisResult

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 20,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr07.csv", row.names=TRUE, quote=FALSE) 
rm(full_results)
rm(plot_list)

#################------------------- Threshold 0.4
#Load target gene list
data = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.4.csv", row.names = 1)
list = data$X0
full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

thisResult1 = data.frame(full_results$results)
thisResult1$MajorCellType = unlist(lapply(strsplit(as.character(thisResult1$CellType),'_'), `[[`, 1))
thisResult1$testList = '0.4 threshold'
FinalResult1 = rbind(FinalResult1,thisResult1)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation1_thr04.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation1_thr04.csv", row.names=TRUE, quote=FALSE) 

rm(full_results)
rm(plot_list)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 2,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length
thisResult = data.frame(full_results$results)
thisResult$MajorCellType = unlist(lapply(strsplit(as.character(thisResult$CellType),'_'), `[[`, 1))
thisResult$testList = '0.4 threshold'
FinalResult = rbind(FinalResult,thisResult)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr04.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 20,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr04.csv", row.names=TRUE, quote=FALSE) 
rm(full_results)
rm(plot_list)

#################------------------- Threshold 0.1
#Load target gene list

data = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.55.csv", row.names = 1)
list = data$X0

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = 50000,                       #Bootstrap repeats set to 10000
                                                annotLevel = 1,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length

thisResult1 = data.frame(full_results$results)
thisResult1$MajorCellType = unlist(lapply(strsplit(as.character(thisResult1$CellType),'_'), `[[`, 1))
thisResult1$testList = '0.55 threshold'
FinalResult1 = rbind(FinalResult1,thisResult1)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation1_thr055.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 7,
  units = c("in"),
  dpi = 300
)

write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation1_thr055.csv", row.names=TRUE, quote=FALSE) 

rm(full_results)
rm(plot_list)

full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                hits = list, 
                                                reps = n,                       #Bootstrap repeats set to 10000
                                                annotLevel = 2,                 #Annotation level where 1= major cell types, 2= higher resolution
                                                geneSizeControl = TRUE)         #Control for GC content and gene length
thisResult = data.frame(full_results$results)
thisResult$MajorCellType = unlist(lapply(strsplit(as.character(thisResult$CellType),'_'), `[[`, 1))
thisResult$testList = '0.55 threshold'
FinalResult = rbind(FinalResult,thisResult)

plot_list = EWCE::ewce_plot(total_res = full_results$results,                     #Write plot data with BH correction p vals
                            mtc_method = "fdr",
                            ctd = ctd)

plot_list$withDendro+ theme(text = element_text(size = 12))+theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))   

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr055.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 20,
  height = 7,
  units = c("in"),
  dpi = 300
)

#write.csv(data.frame(full_results$results), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation2_thr01.csv", row.names=TRUE, quote=FALSE) 
#rm(full_results)
#rm(plot_list)

############-------------------------------- FINAL SAVE and PLOT

write.csv((.packages()), "processed_data/03-LR_network_visualisation/03c-specificity_major_network/loaded_R_packages.R")


write.csv(FinalResult1, "processed_data/03-LR_network_visualisation/03c-specificity_major_network/final_results_table_annot1.csv")
sig1 = FinalResult1[FinalResult1$q>0.05,]

FinalResult1 = read.csv("processed_data/03-LR_network_visualisation/03c-specificity_major_network/final_results_table_annot1.csv", row.names = 1)


FinalResult1 %>% 
  mutate(CellType = factor(CellType,c("Astrocyte", "Microglia","OPC", "Oligodendrocyte","Vascular_cells","GABAergic_LAMP5","GABAergic_PAX6","GABAergic_PVALB","GABAergic_SST","GABAergic_VIP","Glutamatergic_IT","Glutamatergic_L4_IT","Glutamatergic_L5_ET","Glutamatergic_L5_6_IT_Car3","Glutamatergic_L5_6_NP", "Glutamatergic_L6_CT","Glutamatergic_L6b"))) %>%        
  ggplot() +
  geom_point(aes(x = -log10(q), y = CellType,color = MajorCellType, size = fold_change, alpha = q<0.05))+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(cols = vars(testList))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(q)") +
  ylab("Cell type") +
  #coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  guides(fill = guide_legend(title = "Major Cell Type")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.position = "bottom",
        ## This is to plot the two legends in two rows
        legend.box="vertical") 
ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation1_multiple_thr.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 8,
  height = 7,
  units = c("in"),
  dpi = 300
)

sig = FinalResult[FinalResult$q<0.05,]

write.csv(FinalResult, "processed_data/03-LR_network_visualisation/03c-specificity_major_network/annotation1_multiple_thr.csv")

ggplot(FinalResult) +
  geom_point(aes(x = -log10(q), y = CellType, color = MajorCellType, size = fold_change, alpha = q<0.05))+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(rows = vars(testList))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(q)") +
  ylab("Cell type") +
  coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  guides(fill = guide_legend(title = "Major Cell Type")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=270, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.position = "top",
        ## This is to plot the two legends in two rows
        legend.box="vertical") 

ggsave(
  filename="plots/03-LR_network_visualisation/03c-specificity_major_network/annotation2_multiple_thr.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 8,
  height = 7,
  units = c("in"),
  dpi = 300
)
