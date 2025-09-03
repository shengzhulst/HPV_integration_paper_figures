###final_for vis is available at the Supplementary table S6

library(ggpubr)
library(ggplot2)

draw_immuneoncoscore_cor<-function(cell_deconv_final_cor,cell_type,cell_type_name,estimate_use,p_value_use){
    batch_color<-c('UM_FF' = "#F4538A", 'UM_FFPE' = "#FAA300", 'HVC' = "#F5DD61",'TCGA' = "#23f8ff")
    label_used<-paste0('estimate:',estimate_use,'\n','adj p:',p_value_use)
    deconv_vis_immune<-ggplot(cell_deconv_final_cor, aes(full_immunooncoscore_z_score, .data[[cell_type]], color = batch, fill = batch)) +
    geom_point() +
    geom_smooth(method = "lm",se=F) +
    geom_text(aes(x=-3.5,y=max(.data[[cell_type]])-0.03),label=label_used,size=8,color='black') +
    scale_color_manual(values = batch_color, aesthetics = c("color", "fill")) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5),text = element_text(size = 20),axis.text.y=element_text(size=20))+
    labs(title=paste0(cell_type_name),color="Cohort",fill="Cohort")+xlab('Adjusted HGER')+ylab('cell proportion')
    return(deconv_vis_immune)
}

plot_vis_list<-list()

for (i in seq(1:length(cell_type_used))){
    print(i)
    temp_plot<-draw_immuneoncoscore_cor(final_for_vis,cell_type_used[i],cell_type_name[i],estimates_vis[i],p_values_vis[i])
    plot_vis_list[[i]]<-temp_plot
}

plot_vis_list_tertile<-list()
for (i in seq(1:length(cell_type_used))){
    print(i)
    temp_plot<-draw_immunoncoscore_ANOVA_tertile(final_for_vis,cell_type_used[i],cell_type_name[i])
    plot_vis_list_tertile[[i]]<-temp_plot
}

plot_vis_list_tertile[[1]]
plot_vis_list_tertile[[2]]
plot_vis_list_tertile[[3]]
plot_vis_list_tertile[[4]]

library(gridExtra)
library(ggpubr)
png('Figure4C.png',width=17,height=12,units='in',res=300)
ggarrange(plot_vis_list[[1]],plot_vis_list[[2]],plot_vis_list[[3]],plot_vis_list[[4]],nrow = 2,ncol=2,common.legend = TRUE,legend = "right")
dev.off()