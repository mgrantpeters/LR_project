library("biomaRt")
library("ggplot2")
library("grid")
library("dplyr")
set.seed(1234)
setwd("/Users/melis/Documents/GitHub/LR_project/")
data = read.csv("processed_data/06a_splicing_disease_LRs/significance_splicing_disease.csv")

data %>% 
  mutate(Disease = factor(Disease,rev(c("AD","ALS","EssentialTremor","FrontotemporalDementia","LBD","PD","ProgressiveSupranuclearPalsy","AnorexiaNervosa","BipolarDisorder","MajorDepressiveDisorder","NeuroticDisorder","OCD","Schizophrenia","TouretteSyndrome","UnipolarDepression","BrainAneurysm","IntracranialHemorrhage","MigraineDisorder","MigraineWithAura","MS","NarcolepsyCataplexy","Narcolepsy","PartialEpilepsy","RestlessLeg")))) %>%
  ggplot() +
  geom_point(aes(x = -log10(pval), y = Disease,color = Disease, size = Median.transcript.variants, alpha = pval<0.05), show.legend = TRUE)+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(cols = vars(Threshold))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(p)") +
  ylab("Disease") +
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"),
        legend.position = "top")

ggsave(
  filename="plots/06a_splicing_disease_LRs/LR_variants_across_diseases.pdf",
  plot = last_plot(),
  device ="pdf",
  scale = 1,
  width = 10,
  height = 8,
  units = c("in"),
  dpi = 300
)