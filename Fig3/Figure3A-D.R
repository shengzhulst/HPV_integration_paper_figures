library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
library(data.table)
library(jsonlite)


###survival_patitents_vis_no_ND is the table with OS and OS.time and with E1* integration status information 
#all information is avialable at the Supplementary Table 2.1 and Supplementary Table 4.3 and 4.4

visulize_coxph_by_ggplot2<-function(model_for_all,final_title){
summ_all<-summary(model_for_all)$coefficients
summ_all_interval<-summary(model_for_all)$conf.int
number_patients<-summary(model_for_all)$n
final_draw_df<-as.data.frame(cbind(summ_all,summ_all_interval,number_patients))
final_draw_df <- final_draw_df[, !duplicated(as.list(final_draw_df))]
final_draw_df
colnames(final_draw_df)<-c('coef','HR','sd','z','p','negHR','lower','upper','n')
final_draw_df$color_use<-ifelse(final_draw_df$p<0.05,'red','black')
#remove the part inside the `` and the `` for the rownames(final_draw_df)
rownames(final_draw_df)<-gsub("`[^`]*`", "", rownames(final_draw_df))
#for lower or upper, change the inf to 0
final_draw_df$lower[is.infinite(final_draw_df$lower)]<-0
final_draw_df$upper[is.infinite(final_draw_df$upper)]<-0
final_draw_df$name<-factor(rownames(final_draw_df),levels = rev(rownames(final_draw_df)))
final_draw_df$hazard_text<-paste0(round(final_draw_df$HR,digit=3),'\n(',round(final_draw_df$lower,digit=2),'-',round(final_draw_df$upper,digit=2),')')
final_draw_df$HR<-as.numeric(final_draw_df$HR)
max_gap<-max(final_draw_df$upper)-min(final_draw_df$lower)
left_gap<-max_gap/4
right_gap<-max_gap/5
plot_for_hazard_ratio<-ggplot(final_draw_df,aes(y=name,color=color_use))+theme_classic()+
  geom_point(aes(x=HR), shape=15, size=6) +
  geom_errorbar(aes(xmin=lower, xmax=upper),width=0.15)+
  coord_cartesian(ylim=c(1,nrow(final_draw_df)), xlim=c(round(min(final_draw_df$lower))-left_gap-1.5, round(max(final_draw_df$upper))+right_gap+1.5))+
  geom_text(
    aes(x =round(max(final_draw_df$upper))+1+0.5, y = name, label = round(p,digit=4)),
    hjust = 0,size=8)+
  geom_text(
    aes(x = round(min(final_draw_df$lower))-0.5, y = name, label = hazard_text),hjust = 1,size=8)+
   geom_vline(xintercept = 1, linetype="dashed")+ 
  scale_color_manual(values = c('black','red'))+ 
    theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.title.y= element_blank(),text=element_text(size=25),plot.title = element_text(hjust = 0.5),axis.title = element_text(size = 25),legend.position = "none")+
     labs(title = paste0('Multivariable HR(n=',final_draw_df$n,") for ",final_title),x='Hazard Ratio(95%)')
return(plot_for_hazard_ratio)
}

surv_object_full_no_ND <- Surv(time = survival_patitents_vis_no_ND$OS.time , event = survival_patitents_vis_no_ND$OS)

###
#E1* integration status
###
png('Figure3A.png',width = 16, height = 10, units = 'in', res = 300)
survfit2(surv_object_full_no_ND~E1_isoform_integration , data = survival_patitents_vis_no_ND) %>%
      ggsurvfit() +
      labs(
        x = "Months",
        y = "Survival Probability",
        title = 'HPV integration E1* status Kaplan-Meier Plot',
      )  +
      add_confidence_interval() +
      add_risktable(size=9,risktable_stats = c("n.risk"),stats_label = list(n.risk = "Number at Risk",size=12),  theme = theme_risktable_default(axis.text.y.size = 25,
                                                   plot.title.size = 25))+
      add_pvalue("annotation",size=8, pvalue_fun= function(x) style_pvalue(x, digits = 3),caption = "Log-rank {p.value}",x=70)+
      theme_classic()+scale_color_manual(values = c("E1* integration"="#590696","Not E1* integration"= "#FBCB0A", "HPVint(-)"="#7ED7C1"))+
      scale_fill_manual(values = c("E1* integration"="#590696","Not E1* integration"= "#FBCB0A", "HPVint(-)"="#7ED7C1"))+
      theme(text=element_text(size=25),plot.title = element_text(hjust = 0.5),axis.title = element_text(size = 25),
           axis.text = element_text(size = 25),
           legend.text = element_text(size = 25))
dev.off()


survival_final_coxPH_without_HPV_ND_sub1<- survival_patitents_vis_clinical_without_HPV_ND_subset[survival_patitents_vis_clinical_without_HPV_ND_subset$`E1* integration`!='Not E1* integration',]
survival_final_coxPH_without_HPV_ND_sub1$`E1* integration`<-factor(survival_final_coxPH_without_HPV_ND_sub1$`E1* integration`,levels=c('HPVint(-)','E1* integration'))

png('Figure3B.png',width = 14, height = 9, units = 'in', res = 300)
model_for_all <- coxph(Surv(OS.time, OS) ~ `E1* integration`+AJCC8_stage + age + smoker,
                data = survival_final_coxPH_without_HPV_ND_sub1)
model_for_all
visulize_coxph_by_ggplot2(model_for_all,'E1* integration')
dev.off()

#######
##HGER
#######
###recurrence_table_no_ND and recurrence_table_subset_clinical information are also avialable at the Supplementary Table 2.1 and Supplementary Table 4.3 and 4.4
##immune_group_all is the tertile based on HGER
recur_object_full_no_ND<- Surv(time = as.numeric(recurrence_table_no_ND$Recurrence_Time_Months),event = recurrence_table_no_ND$Recurrence_Status)

recur_object_full<- Surv(time = as.numeric(recurrence_table_subset$Recurrence_Time_Months),event = recurrence_table_subset$Recurrence_Status)

png('Figure3C.png',width = 16, height = 10, units = 'in', res = 300)
survfit2(recur_object_full~ immune_group_all, data = recurrence_table_subset) %>%
      ggsurvfit() +
      labs(
        x = "Months",
        y = "Recurrence Free Probability",
        title = 'HGER Tertile FFPE Kaplan-Meier Plot',
      )  +
      add_confidence_interval() +
      add_risktable(size=10,risktable_stats = c("n.risk"),stats_label = list(n.risk = "Number at Risk",size=13),  theme = theme_risktable_default(axis.text.y.size = 25,
                                                   plot.title.size = 25))+
      add_pvalue("annotation",size=9,caption = "Log-rank {p.value}",x=50)+
      theme_classic()+
      theme(text=element_text(size=25),plot.title = element_text(hjust = 0.5),axis.title = element_text(size = 25),
           axis.text = element_text(size = 25),
           legend.text = element_text(size = 25))
dev.off()

png('Figure3D.png',width = 13, height = 9, units = 'in', res = 300)
model_for_all <- coxph(Surv(Recurrence_Time_Months, Recurrence_Status) ~ HGER + age + AJCC8_stage, data = recurrence_table_subset_clinical)
visulize_coxph_by_ggplot2(model_for_all,'HGER')
dev.off()
