library(karyoploteR)
library(data.table)

#######The Figure 1B table is extracted from the Supplmentary Table S1.1
###An example is shown in the directory

temp_for_vis<-fread('Figure1B_table.csv')
#sorte by seqnames and start
temp_for_vis<-temp_for_vis[order(temp_for_vis$seqnames,temp_for_vis$start),]
temp_for_vis
temp_for_vis$col
temp_for_vis<-temp_for_vis %>% group_by(HmGene) %>% arrange(HmGene, y) %>%mutate(y_adjusted = y + cumsum(duplicated(y)) * 0.1) %>%ungroup()

temp_for_vis$color

png('Figure1B.png',width=18,height=12,units = 'in',res = 300)
kp <- plotKaryotype(genome="hg38",chromosomes=unique(as.character(sort(seqnames(temp_for_vis_gr)))),cex=2.5)
kpPoints(kp, chr=temp_for_vis$seqnames, x=temp_for_vis$start, y=temp_for_vis$y_adjusted, pch=temp_for_vis$pch,cex=temp_for_vis$cex, col=temp_for_vis$color)
kpPlotMarkers(kp, data=temp_for_vis_gr, labels=temp_for_vis_gr$HmGene,line.color="white",text.orientation = 'horizontal',max.iter=10000,adjust.label.position=T)
legend(x = "right", pch=c(0,1,2,3),legend = c("UM", "HVC","TCGA",'ScRNA'),cex=1.7)
legend(x = "bottomright", col = c('orange',"red", "green","blue","purple","#00e7fc",'black'), pch=c(0,0,0,0,0,0,0),legend = c('gene_exons',"genes_first_exon", "genes_1to5kb","nearest_gene","genes_introns","enhancer_or_nearest_TSS",'nearby_genes'),cex=1.7)
dev.off()