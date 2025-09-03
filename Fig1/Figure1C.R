library(data.table)
library(dplyr)
library(chipenrich.data)
library(GenomicRanges)
library(data.table)


### z score table
human_final<-read.csv('all_combined_genes_annot_filtered_sporeID_with_z_score.csv',header = T,stringsAsFactors = F)

#############################
#visulize the heatmap for the HPV integration recurrent events
#############################
#all recurrent events(can be find in Supplmentary Table S1)
human_recurrent_final_vis_filtered<-readRDS('Figure1A_human_recurrent_final_vis_filtered.rds')
unique(human_recurrent_final_vis_filtered$SampleID)
#remove ScRNA-seq samples
human_recurrent_final_vis_filtered<-human_recurrent_final_vis_filtered[!grepl("Sc_", human_recurrent_final_vis_filtered$SampleID),]
human_recurrent_final_vis_filtered<-as.data.frame(human_recurrent_final_vis_filtered)
#create sample chr ID to filter
human_recurrent_final_vis_filtered$sample_chr<-paste(human_recurrent_final_vis_filtered$SampleID,human_recurrent_final_vis_filtered$Chr,sep = "_")
human_final$sample_chr<-paste(human_final$SampleID,human_final$Chr,sep = "_")
human_recurrent_final_vis_filtered$left_temp<-human_recurrent_final_vis_filtered$Start-3000000
human_recurrent_final_vis_filtered$right_temp<-human_recurrent_final_vis_filtered$Start+3000000

human_final_recurrent<-merge(human_final,human_recurrent_final_vis_filtered,by.x = c('sample_chr','SampleID','Chr'),by.y = c('sample_chr','SampleID','Chr'))

human_final_recurrent<-unique(human_final_recurrent)

human_final_recurrent<-human_final_recurrent[,c('SampleID','Chr','HmGene.x','annot.x','human_gene_z_score')]
colnames(human_final_recurrent)<-c('SampleID','Chr','HmGene','anno','human_gene_z_score')

#change all other annotation integration if it's not left_one,right_one,left_two,right_two or integration
human_final_recurrent$anno[!human_final_recurrent$anno %in% c('left_one','right_one','left_two','right_two','integration')]<-'integration'

###example of human_final_recurrent

# 459 TCGA-CN-A499-01A  chr3      SOX2-OT integration -0.00509397299922628
# 460 TCGA-CN-A499-01A  chr3      DNAJC19    left_one    0.268946313527423
# 461 TCGA-CN-A499-01A  chr3         FXR1   right_two    0.658045229915567
# 462 TCGA-CN-A499-01A  chr3         SOX2   right_one   -0.111221699483148
# 463 TCGA-CN-A499-01A  chr3    LINC01206    left_two    0.333054183534803


#######you have to manually check to make sure there is no recurrent genes for one sample happened at one chr
#remove some incorrect samples 
human_final_recurrent

#extract the samples integration genes with the highest z_score and keep the genes with the highest z_score

human_final_recurrent$sample_gene<-paste(human_final_recurrent$SampleID,human_final_recurrent$HmGene,sep = "_")


The_integration_part<-human_final_recurrent[human_final_recurrent$anno=='integration',]

The_integration_part<-The_integration_part %>% group_by(SampleID,Chr) %>% 
summarize(human_gene_z_score_new=human_gene_z_score[which.max(abs(human_gene_z_score))],gene_symbol=HmGene[which.max(abs(human_gene_z_score))])

The_right_one<-human_final_recurrent_temp[human_final_recurrent_temp$anno=='right_one',]%>% group_by(SampleID,Chr) %>%
summarize(human_gene_z_score_new=human_gene_z_score[which.max(abs(human_gene_z_score))],gene_symbol=HmGene[which.max(abs(human_gene_z_score))],
sample_gene=sample_gene[which.max(abs(human_gene_z_score))])

human_final_recurrent_temp<-human_final_recurrent_temp[!human_final_recurrent_temp$sample_gene %in% The_right_one$sample_gene,]
dim(The_right_one)

The_left_one<-human_final_recurrent_temp[human_final_recurrent_temp$anno=='left_one',]%>% group_by(SampleID,Chr) %>%
summarize(human_gene_z_score_new=human_gene_z_score[which.max(abs(human_gene_z_score))],gene_symbol=HmGene[which.max(abs(human_gene_z_score))],
sample_gene=sample_gene[which.max(abs(human_gene_z_score))])
dim(The_left_one)

human_final_recurrent_temp<-human_final_recurrent_temp[!human_final_recurrent_temp$sample_gene %in% The_left_one$sample_gene,]
The_left_two<-human_final_recurrent_temp[human_final_recurrent_temp$anno=='left_two',]%>% group_by(SampleID,Chr) %>%
summarize(human_gene_z_score_new=human_gene_z_score[which.max(abs(human_gene_z_score))],gene_symbol=HmGene[which.max(abs(human_gene_z_score))],
sample_gene=sample_gene[which.max(abs(human_gene_z_score))])
dim(The_left_two)

human_final_recurrent_temp<-human_final_recurrent_temp[!human_final_recurrent_temp$sample_gene %in% The_left_two$sample_gene,]
The_right_two<-human_final_recurrent_temp[human_final_recurrent_temp$anno=='right_two',]%>% group_by(SampleID,Chr) %>%
summarize(human_gene_z_score_new=human_gene_z_score[which.max(abs(human_gene_z_score))],gene_symbol=HmGene[which.max(abs(human_gene_z_score))],
sample_gene=sample_gene[which.max(abs(human_gene_z_score))])
dim(The_right_two)

#sort by sampleID chr for all tables
The_integration_part<-The_integration_part[order(The_integration_part$SampleID,The_integration_part$Chr),]
The_left_one<-The_left_one[order(The_left_one$SampleID,The_left_one$Chr),]
The_left_two<-The_left_two[order(The_left_two$SampleID,The_left_two$Chr),]
The_right_one<-The_right_one[order(The_right_one$SampleID,The_right_one$Chr),]
The_right_two<-The_right_two[order(The_right_two$SampleID,The_right_two$Chr),]

head(The_integration_part)

##################################
#complexheatmap part
#########################
library(ComplexHeatmap)
library(chipenrich.data)
library(GenomicRanges)
library(data.table)
library(jsonlite)

z_score_fill<-cbind(The_left_two$human_gene_z_score_new,The_left_one$human_gene_z_score_new,The_integration_part$human_gene_z_score_new,The_right_one$human_gene_z_score_new,The_right_two$human_gene_z_score_new)

gene_matrix<-cbind(The_left_two$gene_symbol,The_left_one$gene_symbol,The_integration_part$gene_symbol,The_right_one$gene_symbol,The_right_two$gene_symbol)
gene_matrix

rownames(z_score_fill)<-paste(The_integration_part$SampleID,The_integration_part$gene_symbol,sep = "_")
colnames(z_score_fill)<-c('left_two','left_one','integrated_assigned','right_one','right_two')
rownames(gene_matrix)<-rownames(z_score_fill)

#sort by z_score_fill by the integrated_assigned and keep index to sort gene_matrix
z_score_fill<-z_score_fill[order(z_score_fill[,3],decreasing = T),]
rownames(z_score_fill)
gene_matrix<-gene_matrix[rownames(z_score_fill),]

#####
library(RColorBrewer)
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
#chr_used<-as.character(seqnames(gr_included))[seq(1,305,by=5)]

z_score_fill<-as.matrix(z_score_fill)
z_score_fill


max_score_index<-which(apply(z_score_fill, 1, function(x) max(abs(x)))>2)
max_score_index
names(max_score_index)

expression_pertubated<-rep('Not_outlier',length(z_score_fill)/5)
expression_pertubated[max_score_index]<-'Outlier'
expression_pertubated
temp_new_df<-data.frame(gene_sample=rownames(z_score_fill),expression=expression_pertubated)

#order by expression
#order gene matrix by temp_new_df
temp_new_df<-temp_new_df[order(temp_new_df$expression,decreasing = T),]
temp_new_df
row_ha = rowAnnotation(expression=temp_new_df$expression,show_annotation_name = FALSE,col = list(expression = c('Outlier' = '#97c028','Not_outlier' = '#e43c96')),
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 14),title_gp=gpar(fontsize=14,fontface = 'bold')))
z_score_fill<-z_score_fill[temp_new_df$gene_sample,]
gene_matrix<-gene_matrix[temp_new_df$gene_sample,]

colnames(z_score_fill)<-c('left second','left first','Assigned gene','right first','right second')

png('/nfs/value/data1/shiting/paper_stuff/HPV_integration_tables_figures/figures/Figure1B.png',width = 9.5,height = 18,units = 'in',res = 300)

Heatmap(z_score_fill,col = colorRamp2(seq(-4, 4, length = 3), c("blue", "#EEEEEE", "red")),
        show_row_names=TRUE, show_column_names=TRUE,cluster_rows = F,cluster_columns = F,row_order=temp_new_df$gene_sample,
        row_names_gp = gpar(fontsize = 13),
        heatmap_legend_param = list(col_fun=col_fun,at = c(-4, 0, 4),
        title = "z_score",labels_gp = gpar(fontsize = 14),title_gp=gpar(fontsize=14,fontface = 'bold')), right_annotation = row_ha, column_names_rot = 45,
        layer_fun = function(j, i, x, y, width, height, fill) {
        v = pindex(z_score_fill, i, j)
        l = abs(v) > 2
        print(l)
        grid.text(sprintf("%s", as.vector(gene_matrix)[l]), x[l], y[l], gp = gpar(fontsize = 13))})

dev.off()
