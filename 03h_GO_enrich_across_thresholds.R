library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
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


results = data.frame(enrich_go$p.adjust,enrich_go$Description,enrich_go$Count)
ordered_results = results[order(results$enrich_go.p.adjust),][1:30,]
ggplot(ordered_results, aes(x=enrich_go.Count, y=reorder(enrich_go.Description, -enrich_go.p.adjust), fill = enrich_go.p.adjust))+
  geom_bar(stat = 'identity')+ 
  theme(text = element_text(size = 25))    

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 13,
  height = 10,
  units = c("in"),
  dpi = 300
)


#########################################################

gene = read.csv("processed_data/03-LR_network_visualisation/louvain_largest_cluster_0.1.csv", row.names = 1)
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


results = data.frame(enrich_go$p.adjust,enrich_go$Description,enrich_go$Count)
ordered_results = results[order(results$enrich_go.p.adjust),][1:30,]
ggplot(ordered_results, aes(x=enrich_go.Count, y=reorder(enrich_go.Description, -enrich_go.p.adjust), fill = enrich_go.p.adjust))+
  geom_bar(stat = 'identity')+ 
  theme(text = element_text(size = 25))    

ggsave(
  filename="plots/03-LR_network_visualisation/03h_GO_enrich_across_thresholds/annotation1_thr07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 13,
  height = 10,
  units = c("in"),
  dpi = 300
)
