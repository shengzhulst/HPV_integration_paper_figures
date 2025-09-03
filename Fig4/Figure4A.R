######################
#read all DEG results
########################
#The DEG results are available at the Supplementary table S5.1-3


library(data.table)
HPV_pos_neg<-fread('edgeR_differential_expression_between_HPV_pos_neg_all_batch.csv')
HPV_IMMUNE_up_down<-fread('edgeR_differential_expression_between_immunooncoscore_high_low_all_batch.csv')
HPV_E1_isoform_all<-fread('edgeR_differential_expression_between_E1_isoform_HPV_all_batch.csv')

#pick out genes have the same direction of the logFC in all batches
head(HPV_E1_isoform)

extract_potentially_differential_genes<-function(HPV_pos_neg,num_sig_needed,logFC_threshold=1){
  HPV_pos_neg<-as.data.frame(HPV_pos_neg)
  all_greater_than_zero <- apply(HPV_pos_neg[,grepl('logFC',colnames(HPV_pos_neg))], 1, function(row) sum(row > logFC_threshold)>=num_sig_needed)
  all_smaller_than_zero <- apply(HPV_pos_neg[,grepl('logFC',colnames(HPV_pos_neg))], 1, function(row) sum(row < (-logFC_threshold))>=num_sig_needed)
  HPV_pos_neg_up<-HPV_pos_neg[all_greater_than_zero,]
  HPV_pos_neg_down<-HPV_pos_neg[all_smaller_than_zero,]
  HPV_pos_neg_up$direction<-'up'
  HPV_pos_neg_down$direction<-'down'
  final_all_potential_genes<-rbind(HPV_pos_neg_up,HPV_pos_neg_down)
  up_neg_number<-c(rep('up',dim(HPV_pos_neg_up)[1]),rep('down',dim(HPV_pos_neg_down)[1]))
  return (list(final_all_potential_genes,up_neg_number))
}
##########################
#create gene lists with same direction of logFC for both heatmap and the gene set enrichment analysis
##############################
logFC_temp_2<-0.58
all_same_direction_significant_genes_E1_isoform_2<-extract_potentially_differential_genes(HPV_E1_isoform,3,logFC_temp_2)[[1]]
#all_same_direction_significant_genes_E6_isoform_2<-extract_potentially_differential_genes(HPV_E6_isoform,3,logFC_temp_2)[[1]]
all_same_direction_significant_genes_HPV_int_pos_neg_2<-extract_potentially_differential_genes(HPV_pos_neg,3,logFC_temp_2)[[1]]
all_same_direction_significant_genes_immunooncoscore_2<-extract_potentially_differential_genes(HPV_IMMUNE_up_down,3,logFC_temp_2)[[1]]
all_same_direction_significant_genes_E1_isoform_all_2<-extract_potentially_differential_genes(HPV_E1_isoform_all,3,logFC_temp_2)[[1]]


#############
#fisher exact test part(extract genes for gene set enrichment analysis)
#############
library(dplyr)

calculate_fisher_and_keep_new_df_with_up_down_gene<-function(HPV_pos_neg_test,genes_kept,FDR_thre){
library(dplyr, help, pos = 2, lib.loc = NULL)
HPV_pos_neg_test<-as.data.frame(HPV_pos_neg_test)
pvalue_columns <- grep('PValue', colnames(HPV_pos_neg_test), value = TRUE)
# Calculate the test_statistic
HPV_pos_neg_test$test_statistic <- -2 * rowSums(log(HPV_pos_neg_test[, pvalue_columns]), na.rm = TRUE)
#HPV_pos_neg_test<-HPV_pos_neg_test %>% mutate(test_statistic = -2 * (log(PValue_HVC) + log(PValue_TCGA) + log(PValue_FF) + log(PValue_FFPE))) 
HPV_pos_neg_test <- HPV_pos_neg_test %>% mutate(p_value_all = 1 - pchisq(test_statistic, df = 8))
HPV_pos_neg_test<-HPV_pos_neg_test %>% mutate(FDR_all = p.adjust(p_value_all, method = "BH"))
HPV_pos_neg_test$Significant<-'not_sig'
HPV_pos_neg_test$Significant[HPV_pos_neg_test$FDR_all<FDR_thre]<-'sig'
HPV_pos_neg_test<-merge(HPV_pos_neg_test,genes_kept,by='genes',all.y=T)
#only keep significant ones
HPV_pos_neg_test_sig<-HPV_pos_neg_test[HPV_pos_neg_test$Significant=='sig',]
print(table(HPV_pos_neg_test_sig$direction))
up_genes<-HPV_pos_neg_test_sig$genes[HPV_pos_neg_test_sig$direction=='up']
down_genes<-HPV_pos_neg_test_sig$genes[HPV_pos_neg_test_sig$direction=='down']
return(list(HPV_pos_neg_test_sig,up_genes,down_genes))
}


FDR_threshold<-0.05
HPV_int_pos_neg_2<-calculate_fisher_and_keep_new_df_with_up_down_gene(HPV_pos_neg,all_same_direction_significant_genes_HPV_int_pos_neg_2,FDR_threshold)
immune_up_down_2<-calculate_fisher_and_keep_new_df_with_up_down_gene(HPV_IMMUNE_up_down,all_same_direction_significant_genes_immunooncoscore_2,FDR_threshold)
#
E1_up_down_2<-calculate_fisher_and_keep_new_df_with_up_down_gene(HPV_E1_isoform,all_same_direction_significant_genes_E1_isoform_2,FDR_threshold)
E1_all_up_down_2<-calculate_fisher_and_keep_new_df_with_up_down_gene(HPV_E1_isoform_all,all_same_direction_significant_genes_E1_isoform_all_2,FDR_threshold)
E1_up_down_2<-E1_all_up_down_2

########################
#pick out the intersect genes
##########################
all_comparisons_list<-list(HPV_int_pos_neg_2,E1_up_down_2,immune_up_down_2)
#reduce and pick out intersect
all_up_genes<-unlist(lapply(all_comparisons_list,function(x) x[[2]]))
all_down_genes<-unlist(lapply(all_comparisons_list,function(x) x[[3]]))

genes_for_all<- c(all_up_genes,all_down_genes)
genes_for_all

all_genes_pos_use_for_vis<-c(HPV_int_pos_neg_2[[2]],HPV_int_pos_neg_2[[3]])
all_genes_immuneonco_for_vis<-c(immune_up_down_2[[3]],immune_up_down_2[[2]])
all_genes_E1_for_vis<-c(E1_up_down_2[[3]],E1_up_down_2[[2]])

genes_HPV_pos_neg_IMMUNE<-intersect(all_genes_pos_use_for_vis,all_genes_immuneonco_for_vis)
genes_HPV_IMMUNE_E1<-intersect(all_genes_immuneonco_for_vis,all_genes_E1_for_vis)
genes_HPV_pos_neg_E1<-intersect(all_genes_pos_use_for_vis,all_genes_E1_for_vis)

genes_HPV_pos_neg_E1
genes_HPV_pos_neg_IMMUNE
genes_HPV_IMMUNE_E1
final_genes_all<-intersect(genes_HPV_pos_neg_IMMUNE,genes_HPV_IMMUNE_E1)
length(intersect(genes_HPV_pos_neg_IMMUNE,genes_HPV_IMMUNE_E1))

#######################
#read the z score separated by batch
#########################
z_score_all_patients<-fread('Figure1B_z_score_all_filtered.csv')
z_score_all_patients<-as.data.frame(z_score_all_patients)
dim(z_score_all_patients)
rownames(z_score_all_patients)<-z_score_all_patients$V1
z_score_all_patients_for_vis<-z_score_all_patients[genes_for_all,]

library(ComplexHeatmap)
#remove the NA row
z_score_all_patients_for_vis<-z_score_all_patients_for_vis[!is.na(z_score_all_patients_for_vis$V1),]
z_score_all_patients_for_vis<-z_score_all_patients_for_vis[,2:dim(z_score_all_patients_for_vis)[2]]

########
#read the patient meta data
########

###############
#finalized the z score
##################
z_score_all_patients_for_vis_1<-z_score_all_patients_for_vis[,patients_meta$samples]

z_score_all_patients_for_vis_1[is.na(z_score_all_patients_for_vis_1)]<-0

##########################
#draw the separated heatmap for each comparisons
##########################
#######
#HPV int pos versus neg
#######
library(ComplexHeatmap)
batch_color<-c('UM_FF' = "#F4538A", 'UM_FFPE' = "#FAA300", 'HVC' = "#F5DD61",'TCGA' = "#23f8ff")
HPVint_status_color<-c('HPVint(+)'="#F24C3D",'HPVint-ND'='#F2BE22','HPVint(-)'='#22A699')
recurrent_or_not_color<-c('recurrent'='#b535d5','non-recurrent'='#0fa395')
outlier_or_not_color<-c('outlier'='#afe721','not_outlier'='#dc4aa6')
patients_meta_HPV_int_pos_neg<-patients_meta[order(patients_meta$batch,patients_meta$HPV_int_status_final),]
#HPV_int_status = anno_simple(as.character(patients_meta$HPV_int_status_final),col=HPVint_status_color)
z_score_all_patients_HPV_pos_neg<-z_score_all_patients_for_vis_1[all_genes_pos_use_for_vis,patients_meta_HPV_int_pos_neg$samples]

column_ha_HPV_pos_neg<-HeatmapAnnotation(cohort=patients_meta_HPV_int_pos_neg$batch,
HPV_int_status=patients_meta_HPV_int_pos_neg$HPV_int_status_final,
col=list(cohort=batch_color,HPV_int_status=HPVint_status_color))



row_ha_HPV_pos_neg =  anno_mark(at = match(final_genes_all,row.names(z_score_all_patients_HPV_pos_neg)), labels = final_genes_all,which='row')

#row_ha_HPV_pos_neg = rowAnnotation(foo = anno_mark(at=c(1,2,3,4), labels =rownames(z_score_all_patients_HPV_pos_neg)[1:4]))
Heatmap_HPV_int_pos<-Heatmap(name='z_score',z_score_all_patients_HPV_pos_neg,show_row_names = TRUE,show_column_names = F,cluster_columns = F,
cluster_rows=F,top_annotation = column_ha_HPV_pos_neg,column_split=patients_meta_HPV_int_pos_neg$batch,
cluster_row_slices = FALSE)+rowAnnotation(foo=row_ha_HPV_pos_neg)
Heatmap_HPV_int_pos

#immune_onco_score
#sort patient_meta to z_score_all_patients_immunooncoscore by the immunooncoscore
patients_meta$full_immunooncoscore_high<-factor(patients_meta$full_immunooncoscore_high,levels = c('low','medium','high'))
patients_meta_immun<-patients_meta[order(patients_meta$batch,patients_meta$full_immunooncoscore_high),]
patients_meta_immun
z_score_all_patients_immunooncoscore<-z_score_all_patients_for_vis_1[all_genes_immuneonco_for_vis,patients_meta_immun$samples]
#z_score_all_patients_immunooncoscore<-z_score_all_patients_immunooncoscore[!is.na(z_score_all_patients_immunooncoscore),]

column_ha_immunooncoscore<-HeatmapAnnotation(cohort=patients_meta_immun$batch,`HGER`=patients_meta_immun$full_immunooncoscore_high,col=list(cohort=batch_color,HGER=c('high'='blue','low'='red','medium'='green')))

row_ha_immunooncoscore = rowAnnotation(foo = anno_mark(at = match(final_genes_all,row.names(z_score_all_patients_immunooncoscore)), labels = final_genes_all))

Heatmap_immunooncoscore<-ComplexHeatmap::Heatmap(name='z_score',z_score_all_patients_immunooncoscore,show_row_names = F,show_column_names = F,cluster_columns = F,
cluster_rows=F,top_annotation = column_ha_immunooncoscore,column_split=patients_meta_immun$batch,right_annotation = row_ha_immunooncoscore)
Heatmap_immunooncoscore

#E1 isoform
patients_meta_E1_isoform<-patients_meta[order(patients_meta$batch,patients_meta$E1_isoform_integration),]
#change the HPVint(-) to Not_E1_isoform_integration
patients_meta_E1_isoform$E1_isoform_integration<-ifelse(patients_meta_E1_isoform$E1_isoform_integration=='HPVint(-)','Not_E1_isoform_integration',patients_meta_E1_isoform$E1_isoform_integration)
patients_meta_E1_isoform<-patients_meta_E1_isoform[patients_meta_E1_isoform$E1_isoform_integration %in% c('E1_isoform_integration','Not_E1_isoform_integration'),]
patients_meta_E1_isoform$E1_isoform_integration<-ifelse(patients_meta_E1_isoform$E1_isoform_integration=='E1_isoform_integration','E1*integration','Not E1*integration')

patients_meta_E1_isoform<-patients_meta_E1_isoform[patients_meta_E1_isoform$HPV_int_status_final %in% c('HPVint(+)','HPVint(-)'),]

table(patients_meta_E1_isoform$E1_isoform_integration,patients_meta_E1_isoform$batch)

E1_isoform_color<-c('E1*integration'="#590696",'Not E1*integration'='#FBCB0A')

z_score_all_patients_E1_isoform<-z_score_all_patients_for_vis_1[all_genes_E1_for_vis,patients_meta_E1_isoform$samples]
head(z_score_all_patients_E1_isoform)

z_score_all_patients_E1_isoform<-z_score_all_patients_E1_isoform[!(row.names(z_score_all_patients_E1_isoform)=='NA'),]
z_score_all_patients_E1_isoform

column_ha_E1_isoform<-HeatmapAnnotation(cohort=patients_meta_E1_isoform$batch,E1star_integration=patients_meta_E1_isoform$E1_isoform_integration,col=list(cohort=batch_color,E1star_integration=E1_isoform_color))
row_ha_E1_isoform = rowAnnotation(foo = anno_mark(at = match(final_genes_all,row.names(z_score_all_patients_E1_isoform)), labels = final_genes_all))

dim(z_score_all_patients_E1_isoform)
Heatmap_E1_isoform<-ComplexHeatmap::Heatmap(name='z_score',z_score_all_patients_E1_isoform,show_row_names = F,show_column_names = F,cluster_columns = F,
cluster_rows=F,top_annotation = column_ha_E1_isoform,column_split=patients_meta_E1_isoform$batch,right_annotation = row_ha_E1_isoform)
Heatmap_E1_isoform

library(gridExtra)
library(ggpubr)
library(patchwork)

p1 = grid.grabExpr(draw(Heatmap_HPV_int_pos))
p2 = grid.grabExpr(draw(Heatmap_immunooncoscore))
p3 = grid.grabExpr(draw(Heatmap_E1_isoform))


png('Figure4A.png',width = 12,height = 23,units = 'in',res = 300)
grid.arrange(p1, p2,p3,
          ncol = 1, nrow = 3)
dev.off()

