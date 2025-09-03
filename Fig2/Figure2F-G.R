library(data.table)
library(ggplot2)
library(unikn)
library(grDevices)
library(ggpubr)
library(patchwork)
dynamic_color<-c(hcl.colors(6, palette = "Spectral"))

########
##### score_E2_E6_E7 and all_samples_immunooncoscore_df information is available at Supplementary Table S2.1, Supplementary Table S4.1, 4.2


score_E2_E6_E7<-fread('all_HPV_samples_int_status_two_judgement_manullay.csv')
colnames(score_E2_E6_E7)
dim(score_E2_E6_E7)
all_samples_immunooncoscore_df<-fread('immune_onco_score_HPV_integration_fusion_all_samples.csv')
dim(all_samples_immunooncoscore_df)

all_score_used_for_vis<-merge(all_samples_immunooncoscore_df,score_E2_E6_E7,by.x="samples",by.y="sample_name",all.x = T)
dim(all_score_used_for_vis)
all_score_used_for_vis<-all_score_used_for_vis[order(all_score_used_for_vis$E2_E6_ratio, decreasing = TRUE),]
head(all_score_used_for_vis)
colnames(all_score_used_for_vis)
fwrite(all_score_used_for_vis,'immune_onco_score_HPV_integration_status_ratio.csv')

#remove all cell_line result
all_score_used_for_vis<-all_score_used_for_vis[all_score_used_for_vis$batch.x!='cell_line',]

hist(all_score_used_for_vis[all_score_used_for_vis$batch.x=='FFPE',]$full_immunooncoscore)
hist(all_score_used_for_vis[all_score_used_for_vis$batch.x=='FFPE',]$partial_immunooncoscore)

all_score_used_for_vis$batch.x[!all_score_used_for_vis$batch.x %in% c('FFPE','single_cell')]<-'FF'
all_score_used_for_vis$whether_FF<-ifelse(all_score_used_for_vis$batch.x=='FF','FF','Not_FF')
all_score_used_for_vis$samples<-factor(all_score_used_for_vis$samples,levels = all_score_used_for_vis$samples)

png('Figure2F.png',width = 18, height = 10, units = 'in', res = 300)
P1<-ggplot(all_score_used_for_vis[all_score_used_for_vis$batch.x=='FF'], aes(samples, E2_E6_ratio, fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "FF", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5))+
        labs(y = "log2(E2/E6)")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free')
P2<-ggplot(all_score_used_for_vis[all_score_used_for_vis$batch.x=='FFPE'], aes(samples, E2_E6_ratio, fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "FFPE", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5))+
        labs(y = "log2(E2/E6)")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free'
P3<-ggplot(all_score_used_for_vis[all_score_used_for_vis$batch.x=='single_cell'], aes(samples, E2_E6_ratio, fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "single cell", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5),legend.position='none')+
        labs(y = "log2(E2/E6)")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free')
layout<-"
AAAAAA
BBBCCC"
P1 + P2 + P3+plot_layout(design = layout,guides = 'collect') 

dev.off()

###
#immuneoncoscore
###
all_score_used_for_vis_immune<-all_score_used_for_vis[order(all_score_used_for_vis$full_immunooncoscore, decreasing = TRUE),]
all_score_used_for_vis_immune$samples<-factor(all_score_used_for_vis_immune$samples,levels = all_score_used_for_vis_immune$samples)


png('Figure2G.png',width = 18, height = 10, units = 'in', res = 300)
P1<-ggplot(all_score_used_for_vis_immune[all_score_used_for_vis_immune$batch.x=='FF'], aes(samples,full_immunooncoscore , fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "FF", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5))+
        labs(y = "HGER")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free')
P2<-ggplot(all_score_used_for_vis_immune[all_score_used_for_vis_immune$batch.x=='FFPE'], aes(samples, full_immunooncoscore, fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "FFPE", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5))+
        labs(y = "HGER")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free'
P3<-ggplot(all_score_used_for_vis_immune[all_score_used_for_vis_immune$batch.x=='single_cell'], aes(samples, full_immunooncoscore, fill =HPV_int_status_final.y )) +
  geom_col() +
  labs(title = "single cell", x = NULL) +
  scale_fill_manual(name='HPVint status',values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D")) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(size=25),
        plot.title = element_text(hjust = 0.5),legend.position='none')+
        labs(y = "HGER")
        #facet_wrap(~batch.x,scales = "free",nrow=3)
        #facet_grid(whether_FF~batch.x,scales = "free",space='free')
layout<-"
AAAAAA
BBBCCC"
P1 + P2 + P3+plot_layout(design = layout,guides = 'collect')

dev.off()