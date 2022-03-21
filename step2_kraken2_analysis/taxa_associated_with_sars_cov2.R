library(tidyverse)
getwd()
<<<<<<< Updated upstream
tib<-read.table(file = "unique_genera2.txt",sep = "\t",header = T)
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
=======

# import the results
tib<-read.table(file = "unique_genera2.txt",sep = "\t",header = T)
tib<-as_tibble(tib)
tib

#filter only those taxa that were signigifacnt to COVID19
tib<-tib%>%filter(treatment_1=="COVID_19")
tib

# make a table with benjamini hochberg adjustments for multiple correction
tib$bh_q_val<-p.adjust(tib$wilcox_p_value, method = "BH", n = length(tib$wilcox_p_value))

#filter out anything with a p adjusted value of > 0.05 
tib<-tib%>%filter(bh_q_val<0.05)
tib

# filter out anything that isnt present at least at a genus level
tib<-tib%>%filter(!is.na(genus))
tib

#Combine the columned treatement 1 and treatment to into COVID vers tmt1 and covidvs tmt2
# and get rid of the HealthyvsCommunity_acquired_pneumonia tab
 tib_tmp<-tib%>%
   group_by(treatment_1,treatment_2,genus,species)%>%
   summarise(across(where(is.numeric),mean))
# 



tib2<-tib_tmp%>%
  pivot_wider(names_from =c(treatment_1,treatment_2),names_sep="vs",values_from=wilcox_p_value)%>%#select(-HealthyvsCommunity_acquired_pneumonia)%>%
  pivot_longer(c(COVID_19vsCommunity_acquired_pneumonia,COVID_19vsUninfected),names_to="comparison",values_to="wilcox_p_value")

tib2
#get rid of any rows without a wilcox_p_value 
#treatments that were significant in 1 comparison ie:(COVID19vstmtX) but not in the other
tib2<-tib2%>%filter(!is.na(wilcox_p_value))
tib2#455 x 16

#arrange the comparisons from the most to least significant
tib2<-tib2%>%arrange(wilcox_p_value)#%>%unnest(wilcox_p_value)
tib2

library(mosaic)
#get rid of the stat significant stuff that only had little bitty changes 
#(or maybe dont)
tib3<-tib2%>%filter(abs(log2_median_ratio)>1) 
tib3

#averaging stuff to look at for fun, dont actually use this
# tib4<-tib3%>%group_by(genus,comparison)%>%summarise(log2_median_ratio=mean(log2_median_ratio),
#                                      median_diff=mean(abs(median_diff)),
#                                      mean_diff=mean(abs(mean_diff)), 
#                                      wilcox_p_value=mean(wilcox_p_value),
#                                      bh_q_val=mean(bh_q_val))%>%filter(!is.na(wilcox_p_value))%>%distinct()
# tib4


#Make a vector that had the most prevalent genera in the list (only greater than 10)
vector<-tally(~genus,tib3,"data.frame")%>%filter(Freq>8)
vector
#filter the most significant genera using the vector
tib5<-tib3%>%
  filter(genus%in%vector$genus)%>%
  filter(!is.na(genus))
tib5
#get rid of the species names and make the genera distinct
# and arrange them by the most significant
tib5<-tib5%>%select(-species)%>%distinct()%>%arrange(genus,wilcox_p_value)
unique(tib5$genus)
unique(tib5$log2_median_ratio)
unique(tib5$mean_diff)
unique(tib5$median_diff)
as.data.frame(tib5)
tib5
write.table(x = tib5,file = "top_taxa_assoc_w_SARS_COV2_infection.tsv",sep = "\t",row.names = F)

>>>>>>> Stashed changes

library(ggpubr)

sphing<-tib%>%filter(genus==c( "Sphingomonas" ))
sphing%>%group_by(genus,treatment_2)%>%summarise(log2_median_ratio=mean(log2_median_ratio),
                                                   median_diff=mean(median_diff),
                                                   mean_diff=mean(mean_diff), 
                                                   wilcox_p_value=mean(wilcox_p_value),
                                                   bh_q_val=mean(bh_q_val))%>%filter(!is.na(wilcox_p_value))%>%distinct()



<<<<<<< Updated upstream
ggerrorplot(data = tib4,x = "log2_median_ratio",y = "genus",color="comparison",numeric.x.axis = T)
write.table(tib4,"top19_genera")
=======
ggerrorplot(data = tib5,x = "log2_median_ratio",y = "genus",color="comparison",numeric.x.axis = T)

>>>>>>> Stashed changes
