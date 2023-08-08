library("EWCE")
library("biomaRt")
library("ggplot2")
library("grid")
library("dplyr")
set.seed(1234)
setwd("/Users/melis/Documents/GitHub/LR_project/")
data = read.csv("processed_data/LR_occurance.csv")
data$padj = p.adjust(data$pval, method="BH")

data %>% 
  mutate(Disease = factor(Disease,rev(c("AD","ALS","EssentialTremor","FrontotemporalDementia","LBD","PD","ProgressiveSupranuclearPalsy","AnorexiaNervosa","BipolarDisorder","MajorDepressiveDisorder","NeuroticDisorder","OCD","Schizophrenia","TouretteSyndrome","UnipolarDepression","BrainAneurysm","IntracranialHemorrhage","MigraineDisorder","MigraineWithAura","MS","NarcolepsyCataplexy","Narcolepsy","PartialEpilepsy","RestlessLeg")))) %>%
  ggplot() +
  geom_point(aes(x = -log10(padj), y = Disease,color = Disease, size = perc.LR, alpha = padj<0.05), show.legend = FALSE)+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(cols = vars(LR))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(p)") +
  ylab("Disease") +
  #coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  #guides(fill = guide_legend(title = "Disease")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
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
  filename="plots/LR_occurance_per_disease_01.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 8,
  units = c("in"),
  dpi = 300
)

#####################################

data = read.csv("processed_data/LR_occurance_0.4.csv")
data$padj = p.adjust(data$pval, method="BH")
data %>% 
  mutate(Disease = factor(Disease,rev(c("AD","ALS","EssentialTremor","FrontotemporalDementia","LBD","PD","ProgressiveSupranuclearPalsy","AnorexiaNervosa","BipolarDisorder","MajorDepressiveDisorder","NeuroticDisorder","OCD","Schizophrenia","TouretteSyndrome","UnipolarDepression","BrainAneurysm","IntracranialHemorrhage","MigraineDisorder","MigraineWithAura","MS","NarcolepsyCataplexy","Narcolepsy","PartialEpilepsy","RestlessLeg")))) %>%
  ggplot() +
  geom_point(aes(x = -log10(padj), y = Disease,color = Disease, size = perc.LR, alpha = padj<0.05), show.legend = FALSE)+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(cols = vars(LR))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(p)") +
  ylab("Disease") +
  #coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  #guides(fill = guide_legend(title = "Disease")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        #legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"))
#legend.title = element_text(colour = "black", size = "10"),
#legend.position = "top",
## This is to plot the two legends in two rows
#legend.box="vertical") 
ggsave(
  filename="plots/LR_occurance_per_disease_04.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 6,
  units = c("in"),
  dpi = 300
)

####################################################

data = read.csv("processed_data/LR_occurance_0.7.csv")
data$padj = p.adjust(data$pval, method="BH")
data %>% 
  mutate(Disease = factor(Disease,rev(c("AD","ALS","PD","BipolarDisorder","Schizophrenia","UnipolarDepression","MigraineDisorder","PartialEpilepsy")))) %>%
  ggplot() +
  geom_point(aes(x = -log10(padj), y = Disease,color = Disease, size = perc.LR, alpha = padj<0.05), show.legend = FALSE)+
  scale_alpha_discrete(range=c(0.4,1))+
  facet_grid(cols = vars(LR))+
  theme(text = element_text(size = 12))+ theme(plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), "in"))+
  geom_vline(xintercept=1.3010, linetype="dashed", color = "darkgray") +
  xlab("-log10(p)") +
  ylab("Disease") +
  #coord_flip() + # this is to make  the graph landscape
  theme_light() +
  viridis::scale_fill_viridis() + ## to make the colour scheme color-blind safe, the parameter option="A" or option="B",etc will change palette
  #guides(fill = guide_legend(title = "Disease")) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "10"),
        axis.text.x = element_text(colour = "black", size = "10", angle=90, hjust=1),
        axis.title = element_text(colour = "black", size = "10"),
        strip.text = element_text(colour = "black", size = "10"), 
        #legend.text = element_text(colour = "black", size = "10"),
        plot.caption = element_text(colour = "black", size = "10"),
        plot.title = element_text(colour = "black", size = "10"))
#legend.title = element_text(colour = "black", size = "10"),
#legend.position = "top",
## This is to plot the two legends in two rows
#legend.box="vertical") 
ggsave(
  filename="plots/LR_occurance_per_disease_07.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 4,
  units = c("in"),
  dpi = 300
)