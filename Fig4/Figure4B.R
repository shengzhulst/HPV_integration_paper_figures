######all_pathway_significant is available at the Supplementary table S5.4-6

library(ggplot2)
library(ggplotify)
#draw ggplot2 dot plots
wrap_label = labeller(comparison = 
    c("HPV_E1_isoform" = "E1*integration",
      "HPV_immune_score" = "HGER",
      "HPV_int_pos_neg" = "HPVintegration\nstatus"))


all_pathway_significant$direction[all_pathway_significant$direction=='up']<-'with'
all_pathway_significant$direction[all_pathway_significant$direction=='down']<-'without'

all_pathway_significant$direction[all_pathway_significant$direction=='with' & all_pathway_significant$comparison=='HPV_immune_score']<-'low'
all_pathway_significant$direction[all_pathway_significant$direction=='without' & all_pathway_significant$comparison=='HPV_immune_score']<-'high'

all_pathway_significant$direction<-factor(all_pathway_significant$direction,levels = c('high','low','without','with'))
all_pathway_significant$direction


#add a row 
all_pathway_significant<-rbind(all_pathway_significant,all_pathway_significant[1,])
all_pathway_significant[1,]$direction<-'without'
all_pathway_significant[1,]$log2_odds_ratio<-NA
all_pathway_significant[1,]$log10_p_adjust<-NA
png('Figure4B.png',width = 20,height =13 ,units = 'in',res = 300)
ggplot(all_pathway_significant)+
  geom_point(aes(x = direction, y = Description, color = log10_p_adjust, size = log2_odds_ratio),stroke=4.5)+
   scale_size(range = c(5, 15)) +
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5),text = element_text(size = 30), axis.title.x = element_blank())+
  labs(x='Direction',
       y='Pathways')+
  coord_cartesian(ylim = c(0, 25))+
  facet_wrap(~comparison,nrow=1,scales = 'free_x',labeller = wrap_label)+
  annotate(
    "segment",
    x = 0.5, xend = 2.5,
    y = 22.5, yend = 22.5,
    arrow = arrow(length = unit(0.6, "cm")),
    color = "red"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 23.5,
    label = "worse clinical \noutcome",
    size = 8,
    hjust = 0.5
  )+
  scale_color_gradient(
    low = "blue", 
    high = "red", 
    name = "-log10(p-adjusted)" # Change the color legend name
  ) +
  scale_size(
    name = "log2(odds ratio)" # Change the size legend name
  )
dev.off()