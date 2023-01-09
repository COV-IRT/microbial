################################################################################
library(tidyverse)
#import the data, index teh values and fix teh canidadtus names
df<-as_tibble(read.table("k:/github/microbial/Additional_file_5_Suppl_Table_6.txt",header = T,sep = "\t"))%>%
  mutate(index=rep(1:length(df$treatment_1)),
         taxon.name=gsub("Candidatus ","Candidatus_",taxon.name))

df2<-df%>%separate(col = taxon.name,remove = F,into = c("genus","trash"),sep = " ",fill = "left")%>%select(-trash)%>%distinct_all()

df2<-df2%>%mutate(genus=gsub("Candidatus_","Candidatus ",genus))
sp<-df2%>%filter(genus=="Sphingomonas")

sp

library(microbiome)
pseq_decontam_no_neg_core<-readRDS("pseq_decontam_no_neg_core.RDS")

sphing<-subset_taxa(pseq_decontam_no_neg_core,genus=="Sphingomonas")
sphing<-subset_samples(sample_sums(sphing)>0,physeq = sphing)
sphing<-psmelt(sphing)
library(ggpubr)
library(ggsci)
my_pal
my_pal2<-c(my_pal,pal_aaas(palette = "default")(10))

names(my_pal2)<-c("COVID19","Community Acquired Pneumonia","Uninfected",unique(sphing$species))
my_pal2
ggboxplot(data = sphing,
          x = "case",
          y = "Abundance",
          color = "case",
          palette = my_pal2,
          title = "Sphingomonas read counts species",
          ylab = "log10 read counts",
          xlab="case cohort",
          add = "jitter",
          add.params = list(size=0.7,color="species"),
          ggtheme = theme_pubr(base_size = 12,legend = "right",x.text.angle = 45))+
  facet_wrap(~species)+
  yscale("log10")

unique(sphing$species)
cap<-df%>%filter(treatment_2=="Community acquired pneumonia")
un<-df%>%filter(treatment_2=="Uninfected")


df%>%filter(taxon.name==contains("Sphingomonas"))
