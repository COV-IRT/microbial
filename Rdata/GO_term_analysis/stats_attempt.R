library(speedyseq)
library(tidyverse)
library(broom)
library(knitr)
a<-psmelt(mol_bac_pseq_prune_deep)%>%
  select(OTU,dmn,Sample,Abundance,case)


case_sum<-psmelt(mol_bac_pseq_prune_deep)%>%
  select(OTU,case,Abundance)%>%
  group_by(OTU,case)%>%summarise(avg=mean(Abundance))%>%
  pivot_wider(id_cols = OTU, names_from = c(case),values_from = avg)
dmm_sum<-psmelt(mol_bac_pseq_prune_deep)%>%
  select(OTU,dmn,Abundance)%>%
  group_by(OTU,dmn)%>%summarise(avg=mean(Abundance))%>%
  pivot_wider(id_cols = OTU, names_from = c(dmn),values_from = avg)


c<-as_tibble(merge(case_sum,dmm_sum))
colnames(c)<-c("OTU",
               "Control_Healthy",
               "Control_Sick",
               "COVID19",
               "dmm1",
               "dmm2",
               "dmm3",
               "dmm4",
               "dmm5")
colnames(c)
library(matrixTests)
library(DESeq2)
dds <- DESeq(dds, minReplicatesForReplace=Inf)



case_sum<-psmelt(mol_bac_pseq_prune_deep)%>%
  select(OTU,case,Abundance)%>%
  group_by(OTU,case)%>%summarise(avg=mean(Abundance))%>%
  pivot_wider(id_cols = case, names_from = OTU,values_from = avg)

krus<-row_kruskalwallis(x = t(case_sum[,2:2387]),case_sum$case)
krus_res<-merge(case_sum,krus)
write.table(krus_res,"krus_case.tsv",sep = "\t")











welch<-row_t_welch(c[,c("COVID19", "Control_Sick", "Control_Healthy")], c[,c("dmm1","dmm2","dmm3","dmm4","dmm5")])

welch<-row_t_welch(c[,c("COVID19", "Control_Sick", "Control_Healthy")], c[,c("dmm1","dmm2","dmm3","dmm4","dmm5")])
d<-merge(c,welch)
d
m1
m2
library(mosaic)
t_test(m1,m2)

case_sum
write.table(c,file = "stats.tsv",sep = "\t")
dim(a)
a
glimpse(a)

stats_df <- a %>% # start with data
  mutate(case = factor(case,levels = c("Control_Healthy", "Control_Sick", "COVID19")),
          dmn= factor(dmn,labels = c("1", "2","3","4","5")))
stats_df
a
df
case_sum
b<-full_join(x = case_sum,y = df,by=Sample)



df<-stats_df%>%pivot_wider(id_cols = c(Sample,case,dmn),names_from=OTU,values_from=Abundance)
write.table(df,file = "stats_df.tsv",sep = "\t")
colnames(stats_df)

# this gives main effects AND interactions
ad_aov <- aov(Abundance ~ OTU*case, data = stats_df)
# this would give ONLY main effects
ad_aov <- aov(Abundance ~ OTU*case*publication*dmn, data = stats_df)


fitted(best)
tmp<-rownames(fitted(best))
test<-as_tibble(fitted(best), row.names=rownames(fitted(best)))
test$bio<-tmp
#test<-pivot_longer(data = test,cols = bio,values_to = c("V1","V2","V3","V4","V5"))
# 
colnames(test)<-c("dmm1","dmm2","dmm3","dmm4","dmm5","bio")
test<-test%>%filter(!is.na(dmm4))
test
glimpse(test)
#ad_aov <- aov(dmm1~ bio,data =test)
# ad_aov
#ad_aov <- aov(bio~ dmm1+dmm2+dmm3+dmm4+dmm5,data =test)
a
rnorm(test[1:5])
normalize <- function(x) { 
  x <- as.matrix(x)
  minAttr=apply(x, 2, min)
  maxAttr=apply(x, 2, max)
  x <- sweep(x, 2, minAttr, FUN="-") 
  x=sweep(x, 2,  maxAttr-minAttr, "/") 
  attr(x, 'normalized:min') = minAttr
  attr(x, 'normalized:max') = maxAttr
  return (x)
} 
test2<-(normalize(a$Abundance))
a
a$norm<-test2

library(mosaic)
test

mat
library(vegan)


logit_mod<-glm(norm~OTU+case,data =a)
bio_aov<-anova(logit_mod)

msummary(logit_mod)
library(broom)
res<-tidy(logit_mod, conf.int =T,exponentiate =T)

tally(~case+dmn,data =a, margins = TRUE)
a
