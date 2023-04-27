library(clusterProfiler)
library(org.Hs.eg.db)
gene = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.7.csv", row.names = 1)
genelist = gene$X0

gene.df <- bitr(genelist, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
sorted = sort(unlist(as.numeric(gene.df$ENTREZID)), decreasing = TRUE)
enrich_go = enrichGO(gene = unlist(sorted),
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  minGSSize = 50,
                  maxGSSize = 300,
                  pvalueCutoff = 0.01,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.01)

barplot(enrich_go, showCategory=20) 

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 4,
  height = 7,
  units = c("in"),
  dpi = 300
)
