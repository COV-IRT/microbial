################################################################################
library(tidyverse)
#get rid of this weird taxa
df<-as_tibble(read.table("k:/github/microbial/Additional_file_5_Suppl_Table_6.txt",header = T,sep = "\t"))#%>%filter(taxon.name!="Chroococcidiopsis thermalis")
df

df<-df%>%mutate(taxon.name=gsub("Candidatus ","Candidatus_",taxon.name))
df$taxon.name
# 
# up<-df%>%filter(log2_median_ratio>0)
# data.frame(up)
# 
# df2<-as_tibble(read.table("k:/github/microbial/step2_kraken2_analysis/statistically_significant_unique_genera.tsv",header = T,sep = "\t"))%>%filter(!is.na(wilcox_p_value))%>%filter(taxon.name!="Chroococcidiopsis thermalis")
# df2
# 
# df2.1<-df2%>%
#   filter(treatment_1=="COVID_19")%>%
#   #filter(abs(log2_median_ratio)>1.0)%>%
#   #filter(!is.na(species))%>%
#   filter(!is.na(genus))%>%
#   filter(!is.na(class))%>%arrange(species)%>%distinct_all()
# 
# df2.1<-df2.1%>%mutate(treatment_2=gsub("Healthy","Uninfected",treatment_2),
#                       treatment_2=gsub("Community_acquired_pneumonia","Community acquired pneumonia",treatment_2),)
# df2.1
# 
# cap<-df%>%filter(treatment_2=="Community acquired pneumonia")%>%arrange(desc(log2_median_ratio))
# 
# un<-df%>%filter(treatment_2=="Uninfected")
# un
# cap
# un%>%filter(log2_median_ratio>0)
# un
# df4<-df2.1%>%filter(!species%in% df$taxon.name)
# df4
# df4<-df2.1%>%filter(!species%in% df$taxon.name)
# df4
# 
# df5<-df%>%filter(!taxon.name%in% df2.1$species)
# df5
# 
# df4
# df$taxon.name
# df
# 

df<-df%>%mutate(taxon.name=gsub("Candidatus ","Candidatus_",taxon.name))
df$taxon.name

cap<-df%>%filter(treatment_2=="Community acquired pneumonia")
un<-df%>%filter(treatment_2=="Uninfected")

cap
un
df2<-df%>%separate(col = taxon.name,into = c("genus","trash"),sep = " ",fill = "left")%>%select(-trash)%>%distinct_all()
df2<-df2%>%mutate(genus=gsub("Candidatus_","Candidatus ",genus))


unique(df2$genus)
df2$genus
sp<-df2%>%filter(genus=="Sphingomonas")
sp

df2

#Average log2 medianfor all species 

df3<-df2%>%group_by(genus,treatment_2,treatment_1)%>%summarise(across(everything(),mean))%>%arrange(genus)%>%ungroup()
df3


gr<-df%>%filter(log2_median_ratio>0)
data.frame(gr)
sp<-df3%>%filter(genus=="Sphingomonas")
sp

covid_un<-df3%>%filter(treatment_2=="Uninfected")%>%
  filter(abs(log2_median_ratio)>2)%>%
  filter(abs(median_diff)>0.05)%>%
  filter(abs(mean_diff)>0.1)%>%
  filter(bh_q_value<0.01)
data.frame(covid_un)

covid_cap<-df3%>%filter(treatment_2=="Community acquired pneumonia")%>%
  filter(abs(log2_median_ratio)>2)%>%
  filter(abs(median_diff)>0.05)%>%
  filter(abs(mean_diff)>0.1)%>%
  filter(bh_q_value<0.01)
data.frame(covid_cap)
covid_cap<-df3%>%filter(treatment_2=="Community acquired pneumonia")
  
covid_cap_slim<-covid_cap%>%filter(genus%in%covid_un$genus)
covid_un_slim<-covid_un%>%filter(genus%in%covid_cap$genus)

both<-full_join(covid_cap_slim,covid_un_slim)



covid_cap_unslim<-covid_cap%>%filter(!genus%in%covid_un$genus)
covid_un_unslim<-covid_un%>%filter(!genus%in%covid_cap$genus)

both_unslim<-full_join(covid_cap_unslim,covid_un_unslim)

data.frame(both_unslim)



covid_cap
covid_un
both
data.frame(covid_un)

data.frame(un)%>%
  filter(bh_q_value<=favstats(un$bh_q_value)$median)%>%
  filter(log2_median_ratio>=favstats(un$log2_median_ratio)$median)%>%arrange(desc(log2_median_ratio))


data.frame(un)%>%slice_min(order_by = log2_median_ratio,n = 5)
data.frame(un)%>%slice_min(order_by = log2_median_ratio,n = 10)

cap<-covid_cap%>%mutate(q.value=p.adjust(p = wilcox_p_value,n = length(cap$wilcox_p_value),method = "hochberg"))%>%filter(q.value<0.05)



covid_cap_slim<-covid_cap%>%filter(genus%in%covid_un$genus)
covid_un_slim<-covid_un%>%filter(genus%in%covid_cap$genus)
df3
df4_min<-df3%>%group_by(treatment_2)%>%slice_min(order_by = log2_median_ratio,n = 1)
df4_min
df4_max<-df3%>%group_by(treatment_2)%>%slice_max(order_by = log2_median_ratio,n = 1)
df4_max
df4<-full_join(df4_min,df4_max)
df4
data.frame(covid_un_slim)
data.framecovid_un_slim

ggscatterhist(df3,x = "median_diff",y = "genus",color = "treatment_2",fill = "treatment_2")+scale_color_aaas()+scale_fill_aaas()
plot(x = df3$median_diff,
     y = rep.int(1,times = length(df3$median_diff)))
plot(df3$median_diff)

# df3<-df2%>%group_by(treatment_1,treatment_2,genus)%>%slice_max(order_by = abs(log2_median_ratio),n = 1)%>%distinct_all()
unique(df3$genus)


library(ggpubr)
library(mosaic)
abs(favstats(df$median_diff)$median)
unique(df$log2_median_ratio)

gghistogram(data = df,x = "mean_diff",bins = 50)
hist(df$mean_diff)

range(notinf$log2_median_ratio)
tail(notinf)
covid_cap<-df%>%filter(treatment_2=="Community acquired pneumonia")
range(data.frame(covid_cap$log2_median_ratio))



covid_un<-df3%>%filter(treatment_2=="Uninfected")
covid_cap_slim<-covid_cap%>%filter(genus%in%covid_un$genus)
covid_un_slim<-covid_cap%>%filter(genus%in%covid_cap$genus)


res<-full_join(covid_cap_slim,covid_un_slim)%>%ungroup%>%arrange(desc(log2_median_ratio))
unique(res$genus)
inf<-res%>%filter(log2_median_ratio=="Inf")
notinf<-res%>%filter(log2_median_ratio!="Inf")
range(notinf$log2_median_ratio)
inf
notinf
