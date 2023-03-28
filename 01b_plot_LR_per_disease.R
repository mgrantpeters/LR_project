library("EWCE")
library("biomaRt")
library("ggplot2")
library("grid")
library("dplyr")
set.seed(1234)
setwd("/Users/melis/Documents/GitHub/LR_project/")
data = read.csv("processed_data/LR_occurance.csv")

data %>% 
  ggplot() +
  geom_point(aes(x = -log10(pval), y = Disease,color = Disease, size = perc.LR, alpha = pval<0.05), show.legend = FALSE)+
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
  filename="plots/LR_occurance_per_disease.png",
  plot = last_plot(),
  device ="png",
  scale = 1,
  width = 7,
  height = 8,
  units = c("in"),
  dpi = 300
)