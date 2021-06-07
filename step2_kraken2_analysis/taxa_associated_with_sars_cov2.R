library(tidyverse)
getwd()
tib<-read.table(file = "unique_genera.tsv",sep = "\t",header = T)
tib<-as_tibble(tib)
tib
#filter only those taxa that were signigifacnt to COVID19
tib<-tib%>%filter(treatment_1=="COVID_19")
tib
#tib<-tib%>%pivot_longer(cols = c(treatment_1,treatment_2),names_to="comparison",values_to="treatment")
tib$bh_q_val<-p.adjust(tib$wilcox_p_value, method = "BH", n = length(tib$wilcox_p_value))
tib<-tib%>%filter(bh_q_val<0.05)
tib

tib<-tib%>%filter(!is.na(genus))

#tib<-tib%>%filter(!is.na(genus))
#tib<-tib%>%select(-c(domain,phylum,class,order,family,genus))
tib
#tib<-tib%>%filter(duplicated(tib$genus))%>%arrange(genus)
#tib<-tib%>%filter(!is.na(log2_median_ratio))
tib

#tib<-tib%>%filter(abs(mean_diff)>0.01)

tib2<-tib%>%pivot_wider(names_from =c(treatment_1,treatment_2),names_sep="vs",values_from=wilcox_p_value)
tib2
tib2<-tib2%>%pivot_longer(c(COVID_19vsCommunity_acquired_pneumonia,COVID_19vsHealthy),names_to="comparison",values_to="wilcox_p_value")
tib2<-tib2%>%filter(!is.na(wilcox_p_value))
tib2
tib2<-tib2%>%arrange(desc(wilcox_p_value))%>%unnest(wilcox_p_value)
tib2

library(mosaic)
tib3<-tib2%>%filter(abs(log2_median_ratio)>1)
tib3
#tib3<-tib3%>%group_by(genus,comparison)%>%summarise(log2_median_ratio=mean(log2_median_ratio),
                                     median_diff=mean(abs(median_diff)),
                                     mean_diff=mean(abs(mean_diff)), 
                                     wilcox_p_value=mean(wilcox_p_value),
                                     bh_q_val=mean(bh_q_val))%>%filter(!is.na(wilcox_p_value))%>%distinct()
tib3
vector<-tally(~genus,tib3,"data.frame")%>%filter(Freq>10)

tib4<-tib3%>%filter(genus%in%vector$genus)
tib4<-tib4%>%select(-c(taxon_id,otu_id))%>%filter(!is.na(genus))
tib4<-tib4%>%select(-species)%>%distinct()
tib4
tib4%>%arrange(genus,wilcox_p_value)
unique(tib4$genus)
unique(tib4$log2_median_ratio)
unique(tib4$mean_diff)
unique(tib4$median_diff)
as.data.frame(tib4)

library(ggpubr)

sphing<-tib%>%filter(genus==c( "Sphingomonas" ))
sphing%>%group_by(genus,treatment_2)%>%summarise(log2_median_ratio=mean(log2_median_ratio),
                                                   median_diff=mean(median_diff),
                                                   mean_diff=mean(mean_diff), 
                                                   wilcox_p_value=mean(wilcox_p_value),
                                                   bh_q_val=mean(bh_q_val))%>%filter(!is.na(wilcox_p_value))%>%distinct()



ggerrorplot(data = tib4,x = "log2_median_ratio",y = "genus",color="comparison",numeric.x.axis = T)
write.table(tib4,"top19_genera")
