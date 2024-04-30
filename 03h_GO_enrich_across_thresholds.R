library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
setwd("/Users/melis/Documents/GitHub/LR_project")
gene = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.7.csv", row.names = 1)
genelist = gene$X0
universe = read.csv("processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/universe_GOenrich_brain_proteins.csv")
universe = universe$gene

gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
universe.df <- bitr(universe, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
sorted = sort(unlist(as.numeric(gene.df$ENTREZID)), decreasing = TRUE)
enrich_go = enrichGO(gene = unlist(sorted),
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     minGSSize = 50,
                     maxGSSize = 300,
                     pvalueCutoff = 0.01,
                     pAdjustMethod = "fdr",
                     universe = universe.df$ENTREZID, 
                     qvalueCutoff = 0.01)

enrich_go = clusterProfiler::simplify(enrich_go,
                                      cutoff = 0.7,
                                      by = "p.adjust")
results = data.frame(enrich_go$p.adjust,enrich_go$Description,enrich_go$Count)
ordered_results = results[order(results$enrich_go.p.adjust),][1:30,]
#ggplot(ordered_results, aes(x=enrich_go.Count, y=reorder(enrich_go.Description, -enrich_go.p.adjust), fill = enrich_go.p.adjust))+
#  geom_bar(stat = 'identity')+ 
#  theme(text = element_text(size = 25))   

# Initialize empty vectors for numerators and denominators
numerators <- numeric(length(enrich_go$GeneRatio))
denominators <- numeric(length(enrich_go$GeneRatio))

# Extract numerators and denominators
for (i in 1:length(enrich_go$GeneRatio)) {
  parts <- strsplit(enrich_go$GeneRatio[i], "/")[[1]]
  numerators[i] <- as.numeric(parts[1])
  denominators[i] <- as.numeric(parts[2])
}
# Calculate the ratios
ratios <- numerators / denominators
enrich_go$Ratios = list(ratios)

results = data.frame(Pathway = enrich_go$Description, GeneRatio = enrich_go$GeneRatio, Counts = enrich_go$Count, Ratios = ratios, P.adjust = enrich_go$p.adjust)
results = results[order(results$P.adjust), ]

top = head(results, 15)
top = top[order(top$Ratios), ]
top$Pathway <- factor(top$Pathway, levels = top$Pathway)

ggplot(top) +
  geom_point(aes(x = Ratios, y = Pathway, color = P.adjust, size = Counts ))+
  scale_color_continuous(low = "cyan3", high = "black")+
  xlab("Gene Ratio") +
  ylab("Gene ontology biological process pathway") +
  theme_light() +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=0, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.box="horizontal",
        legend.direction="vertical",
        legend.position="bottom") 
write.csv(top, file = 'processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.csv')
ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.pdf",
  plot = last_plot(),
  device ="pdf",
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300
)

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 6,
  height = 5,
  units = c("in"),
  dpi = 300
)
#########################################################

gene = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.55.csv", row.names = 1)
genelist = gene$X0
universe = read.csv("processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/universe_GOenrich_brain_proteins.csv")
universe = universe$gene

gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
universe.df <- bitr(universe, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
sorted = sort(unlist(as.numeric(gene.df$ENTREZID)), decreasing = TRUE)
enrich_go = enrichGO(gene = unlist(sorted),
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     minGSSize = 50,
                     maxGSSize = 300,
                     pvalueCutoff = 0.01,
                     pAdjustMethod = "fdr",
                     universe = universe.df$ENTREZID,
                     qvalueCutoff = 0.01)

enrich_go = clusterProfiler::simplify(enrich_go,
                                      cutoff = 0.7,
                                      by = "p.adjust")
results = data.frame(enrich_go$p.adjust,enrich_go$Description,enrich_go$Count)
ordered_results = results[order(results$enrich_go.p.adjust),][1:30,]
#ggplot(ordered_results, aes(x=enrich_go.Count, y=reorder(enrich_go.Description, -enrich_go.p.adjust), fill = enrich_go.p.adjust))+
#  geom_bar(stat = 'identity')+ 
#  theme(text = element_text(size = 25))   

# Initialize empty vectors for numerators and denominators
numerators <- numeric(length(enrich_go$GeneRatio))
denominators <- numeric(length(enrich_go$GeneRatio))

# Extract numerators and denominators
for (i in 1:length(enrich_go$GeneRatio)) {
  parts <- strsplit(enrich_go$GeneRatio[i], "/")[[1]]
  numerators[i] <- as.numeric(parts[1])
  denominators[i] <- as.numeric(parts[2])
}
# Calculate the ratios
ratios <- numerators / denominators
enrich_go$Ratios = list(ratios)

results = data.frame(Pathway = enrich_go$Description, GeneRatio = enrich_go$GeneRatio, Counts = enrich_go$Count, Ratios = ratios, P.adjust = enrich_go$p.adjust)
results = results[order(results$P.adjust), ]

top = head(results, 15)
top = top[order(top$Ratios), ]
top$Pathway <- factor(top$Pathway, levels = top$Pathway)

ggplot(top) +
  geom_point(aes(x = Ratios, y = Pathway, color = P.adjust, size = Counts ))+
  scale_color_continuous(low = "cyan3", high = "black")+
  xlab("Gene Ratio") +
  ylab("Gene ontology biological process pathway") +
  theme_light() +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=0, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.box="vertical") 

write.csv(top, file = 'processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr055.csv')
ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr055.pdf",
  plot = last_plot(),
  device ="pdf",
  scale = 1,
  width = 7,
  height = 5,
  units = c("in"),
  dpi = 300
)

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr055.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 5,
  units = c("in"),
  dpi = 300
)

#########################################################

gene = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.4.csv", row.names = 1)
genelist = gene$X0
universe = read.csv("processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/universe_GOenrich_brain_proteins.csv")
universe = universe$gene

gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
universe.df <- bitr(universe, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
sorted = sort(unlist(as.numeric(gene.df$ENTREZID)), decreasing = TRUE)
enrich_go = enrichGO(gene = unlist(sorted),
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     minGSSize = 50,
                     maxGSSize = 300,
                     pvalueCutoff = 0.01,
                     pAdjustMethod = "fdr",
                     universe = universe.df$ENTREZID,
                     qvalueCutoff = 0.01)

enrich_go = clusterProfiler::simplify(enrich_go,
                                      cutoff = 0.7,
                                      by = "p.adjust")
results = data.frame(enrich_go$p.adjust,enrich_go$Description,enrich_go$Count)
ordered_results = results[order(results$enrich_go.p.adjust),][1:30,]
#ggplot(ordered_results, aes(x=enrich_go.Count, y=reorder(enrich_go.Description, -enrich_go.p.adjust), fill = enrich_go.p.adjust))+
#  geom_bar(stat = 'identity')+ 
#  theme(text = element_text(size = 25))   

# Initialize empty vectors for numerators and denominators
numerators <- numeric(length(enrich_go$GeneRatio))
denominators <- numeric(length(enrich_go$GeneRatio))

# Extract numerators and denominators
for (i in 1:length(enrich_go$GeneRatio)) {
  parts <- strsplit(enrich_go$GeneRatio[i], "/")[[1]]
  numerators[i] <- as.numeric(parts[1])
  denominators[i] <- as.numeric(parts[2])
}
# Calculate the ratios
ratios <- numerators / denominators
enrich_go$Ratios = list(ratios)

results = data.frame(Pathway = enrich_go$Description, GeneRatio = enrich_go$GeneRatio, Counts = enrich_go$Count, Ratios = ratios, P.adjust = enrich_go$p.adjust)
results = results[order(results$P.adjust), ]

top = head(results, 15)
top = top[order(top$Ratios), ]
top$Pathway <- factor(top$Pathway, levels = top$Pathway)
ggplot(top) +
  geom_point(aes(x = Ratios, y = Pathway, color = P.adjust, size = Counts ))+
  scale_color_continuous(low = "cyan3", high = "black")+
  xlab("Gene Ratio") +
  ylab("Gene ontology biological process pathway") +
  theme_light() +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=0, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.title = element_text(colour = "black", size = "10"),
        legend.box="horizontal",
        legend.direction="vertical",
        legend.position="bottom") 

write.csv(top, file = 'processed_data/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr04.csv')

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr04.pdf",
  plot = last_plot(),
  device ="pdf",
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300
)

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr04.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 5,
  height = 5,
  units = c("in"),
  dpi = 300
)