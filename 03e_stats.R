setwd("~/GitHub/LR_project")
res = c('0.4','0.55','0.7')
cluster = c('L1 ~ SpD1','L1 ~ SpD2','L2 ~ SpD3','L5 ~ SpD4','L3 ~ SpD5','WM ~ SpD6','L6 ~ SpD7','L4 ~ SpD8','WM ~ SpD9')
for (resolution in res) {
  pval <- numeric(9)
  for (i in 1:9){
    data = read.csv(sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/scores_for_stats_%s_%s.csv", resolution, cluster[i]))
    result = wilcox.test(data$genes_of_interest_final, data$background_final, paired=TRUE)
    print(result)
    pval[i] = result$p.value
  }
  df <- data.frame(Column1 =cluster, Column2 =pval)
  write.csv(df, sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/stats_results_wilcoxonpaired_%s_%s.csv", resolution, cluster[i]))
}

data = read.csv(sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/res_0.4_scores_genes_interest.csv"))
kruskal.test(score_final ~ bayesSpace_harmony_9, data = data)
res = dunnTest(score_final ~ bayesSpace_harmony_9, data=data, method="bh")
write.csv(res$res, sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/dunn_test_0.4.csv"))


data = read.csv(sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/res_0.55_scores_genes_interest.csv"))
kruskal.test(score_final ~ bayesSpace_harmony_9, data = data)
res = dunnTest(score_final ~ bayesSpace_harmony_9, data=data, method="bh")
write.csv(res$res, sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/dunn_test_0.55.csv"))


data = read.csv(sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/res_0.7_scores_genes_interest.csv"))
kruskal.test(score_final ~ bayesSpace_harmony_9, data = data)
res = dunnTest(score_final ~ bayesSpace_harmony_9, data=data, method="bh")
write.csv(res$res, sprintf("processed_data/03-LR_network_visualisation/03e_network_spatial_analysis/dunn_test_0.7.csv"))
