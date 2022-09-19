library(tidyverse)
df<-as_tibble(read.table("k:/github/microbial/Additional_file_5_Suppl_Table_6.txt",header = T,sep = "\t"))

df<-df%>%mutate(taxon.name=gsub("Candidatus ","",taxon.name))
df_slim<-df%>%filter(abs(log2_median_ratio)>1.0)
df_slim
unique(df_slim$taxon.name)

df2<-df%>%separate(col = taxon.name,into = c("genus","trash"),fill = "left")%>%select(-trash)%>%distinct_all()
length(unique(df2$genus))

sp<-df2%>%filter(genus=="Sphingomonas")
sp


df3<-df2%>%group_by(treatment_1,treatment_2,genus)%>%slice_max(order_by = abs(log2_median_ratio),n = 1)%>%distinct_all()
unique(df3$genus)
inf

library(ggpubr)
unique(df$median_diff)
unique(df$log2_median_ratio)

gghistogram(data = df,x = "mean_diff",bins = 50)
hist(df$mean_diff)

range(notinf$log2_median_ratio)
tail(notinf)
covid_cap<-df3%>%filter(treatment_2=="Community acquired pneumonia")
covid_un<-df3%>%filter(treatment_2=="Uninfected")
covid_cap_slim<-covid_cap%>%filter(genus%in%covid_un$genus)
covid_un_slim<-covid_cap%>%filter(genus%in%covid_cap$genus)
res<-full_join(covid_cap_slim,covid_un_slim)%>%ungroup%>%arrange(desc(log2_median_ratio))
unique(res$genus)
inf<-res%>%filter(log2_median_ratio=="Inf")
notinf<-res%>%filter(log2_median_ratio!="Inf")
range(notinf$log2_median_ratio)
