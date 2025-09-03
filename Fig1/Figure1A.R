#proportion figures
#read the meta data for draw the figure
library(data.table)
library(ggplot2)
##########the meta_data HPV int status final is available in Supplementary Table S2.1
meta_data<-fread('all_HPV_samples_int_status_two_judgement_manullay.csv')

head(meta_data)
table(meta_data$HPV_int_status_final)
meta_data$value<-1
table(meta_data$HPV_int_status_final)
meta_data$HPV_int_status_final<-factor(meta_data$HPV_int_status_final,levels = c('HPVint(-)', 'HPVint-ND','HPVint(+)'))
#calculate the proportion
library(dplyr)
meta_data <- meta_data %>%
  group_by(batch, HPV_int_status_final) %>%
  summarise(value = sum(value)) %>%
  mutate(total = sum(value)) %>%
  ungroup() %>%
  mutate(proportion = value / total)
# meta_data
totals <- meta_data %>%
  group_by(batch) %>%
  summarise(total_batch = sum(value))
totals
#change the FF to UM_FF
meta_data

meta_data$batch[meta_data$batch=='FF']<-'UM_FF'
meta_data$batch[meta_data$batch=='FFPE']<-'UM_FFPE'
####visulization
png('Figure1A.png',width = 3500,height = 3500,units = 'px',res = 300)
ggplot(meta_data, aes(x = batch, y = value, fill = HPV_int_status_final)) +
  geom_col(position = "fill", stat = "identity") + xlab('cohort')+ylab('proportion')+
  geom_text(data=meta_data,aes(label = scales::percent(proportion,accuracy=1L), y = proportion), 
            position = position_fill(vjust = 0.5), size = 7)+ 
            theme_minimal()+
  scale_fill_manual(values = c("HPVint(-)" = "#22A699", "HPVint-ND" = "#F2BE22", "HPVint(+)" = "#F24C3D"))+
  geom_text(data=meta_data,aes(x = batch, label = total, y = 1.03,color='Number of samples'), vjust = 0, size = 7,fontface = "bold")+
  theme(text = element_text(size = 30),plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 45))+
  guides(fill=guide_legend(title="HPV integration status"),color=guide_legend(title="Number of samples"))
dev.off()










