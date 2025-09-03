library(dplyr)
library(reshape2)
library(grid)
library(scales)
library(ggrepel)
library(ggplot2)
library(patchwork)
#### the human final information is available from the Supplementary Table S1 

human_final<-read.csv('all_combined_genes_annot_filtered.csv',header = T,stringsAsFactors = F)

double_results<-human_final %>% group_by(HPVStart) %>% summarise(sample_chr_start_number=n_distinct(sample_chr_start),sample_number=n_distinct(SampleID))

double_results_vis<-melt(double_results,id.vars='HPVStart')
double_results_vis<-double_results_vis[double_results_vis$variable=='sample_number',]

splice_sites<-as.data.frame(c(227,408,526,881,1302,2708,3358,5638))
colnames(splice_sites)<-'sites'

# Input genomic data
hpv16_genes <- data.frame(
  HPV_gene = c("E6", "E7", "E1", "E2", "E5", "L1", "L2"),
  start = c(82,  561, 864, 2754,  3871, 5638, 4235),
  end = c(559,  858, 2813, 3852,  4100, 7149, 5657)
)

# Calculate the length of each gene
hpv16_genes <- hpv16_genes %>%
  mutate(length = end - start + 1)

# Calculate ratios
hpv16_genes <- hpv16_genes %>%
  mutate(ratio = length / sum(length))

# Assign colors manually or generate automatically
hpv16_genes <- hpv16_genes %>%
  mutate(color = alpha(c(
    "orange", "yellow", "red", "blue", "green", "purple", "pink"
  ), 0.3))  # Adjust transparency with alpha

# Prepare for ggplot2
hpv16_genes <- hpv16_genes %>%
  mutate(start_ratio = cumsum(c(0, head(ratio, -1))),
         end_ratio = cumsum(ratio))
# Create the base layer for the bottom bar
bottom_bar<-ggplot(hpv16_genes) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.05, ymax = 0,
                fill = color),
            color = "black") +  # Optional borders
  geom_text(aes(x = (start + end) / 2, y = -0.025, label = HPV_gene),
            size = 9, color = "black") +
  scale_x_continuous( expand = c(0,0) , limits = c(0,7843) )+
  scale_fill_identity() +  # Use predefined fill colors
  theme_void()  # Remove grid and background
bottom_bar

main_plot<-ggplot(double_results_vis) + geom_hline(aes(yintercept=1.5849, color="3 samples"), size=1)+
geom_text_repel(data=splice_sites, aes(x=sites, y=8, label=paste0('SJ:',"\n",sites)), vjust=0.8, color="#4793AF",size=10)+ 
geom_text(aes(x=4181, y=8, label=paste0('TATA:',"\n","4181")), vjust=0.8, color="#4793AF",size=10)+
geom_segment(data=splice_sites,aes(x=sites,y=0,yend=6.5,color="HPV splice junction sites"), size=1)+
geom_segment(aes(x = 4181,y=0,yend=6),color="#4793AF",size=1)+
scale_x_continuous( expand = c(0,0) , limits = c(0,7843) )+
scale_y_continuous( expand = c(0,0), limits = c(0,8.2) )+
geom_point(aes(x = HPVStart, y = log2(value)),size=6)+
theme_minimal()+theme(text=element_text(size = 35),legend.title = element_blank(),legend.position='top')+labs(x='HPV',y='log2(number of samples)')+scale_color_manual(values=c('3 samples'='#DD5746', 'HPV splice junction sites'='#4793AF'))
final_plot <- main_plot / bottom_bar + plot_layout(heights = c(5, 1)) 
final_plot

png('Figure2B.png',width=27,height=10,units='in',res=300)
final_plot
dev.off()   