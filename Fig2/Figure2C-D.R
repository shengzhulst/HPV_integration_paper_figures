#the all_score_used_for_vis_immune information is available at the Supplementary Table S4.1 and S4.2 and Supplementary Table S2.1


my_comparisons <- list( c("HPVint(-)", "E1* integration"), c("Not E1* integration", "HPVint(-)"), c("E1* integration", "Not E1* integration") )
png('Figure2C.png',width = 35, height = 15, units = 'in', res = 300)
ggplot(all_score_used_for_vis_immune,aes(x=E1_isoform_integration,y=E2_E6_ratio,fill=E1_isoform_integration))+
  geom_violin()+
  stat_compare_means(comparisons=my_comparisons,size=13)+
  geom_boxplot(width=0.1)+
  theme_classic()+
  scale_fill_manual(values=c("#590696", "#FBCB0A", "#7ED7C1"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.line.x = element_blank(),
        text = element_text(size = 45,family='Arial'),legend.title=element_blank())+
  labs(y='log2(E2/E6)')+
  facet_wrap(~batch.y,nrow=1)
dev.off()

my_comparisons <- list( c("HPVint(-)", "E6* integration"), c("Not E6* integration", "HPVint(-)"), c("E6* integration", "Not E6* integration") )
png('Figure2D.png',width = 35, height = 15, units = 'in', res = 300)
ggplot(all_score_used_for_vis_immune,aes(x=E6_isoform_integration,y=E2_E6_ratio,fill=E6_isoform_integration))+
  geom_violin()+
  stat_compare_means(comparisons=my_comparisons,size=13)+
  geom_boxplot(width=0.1)+
  theme_classic()+
  scale_fill_manual(values=c("#B06161", "#F0DBAF", "#7ED7C1"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), axis.line.x = element_blank(),
        text = element_text(size = 45),legend.title=element_blank())+
  labs(y='log2(E2/E6)')+
  facet_wrap(~batch.y,nrow=1)
dev.off()
