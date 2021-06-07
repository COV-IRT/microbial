setwd('/home/jovyan/work/microbial/GO_term_analysis/16_DEC_2020_dmm_troubleshooting_Lee')

library(tidyverse)
library(phyloseq)
setwd("./")

raw<-as_tibble(read.table("../Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)

df$depth<-as.character(df$depth)

term<-df%>%filter(type!="NA")%>%filter(type%in%c("bac","arc","vir"))%>%group_by(GO_term,namespace,depth,name)%>%
  summarise(across(.cols = where(is.numeric), sum))

term_tax<-term%>%select(GO_term,namespace,depth,name)
term_tax<-data.frame(term_tax, row.names=1)
term_counts<-data.frame(term[5:172], row.names = term$GO_term)

term_counts_phy <- otu_table(term_counts, taxa_are_rows=TRUE)
term_tax_phy <- tax_table(as.matrix(term_tax), errorIfNULL=TRUE)

term_sam<-as.data.frame(read.table("../Combined_BALF_GO_Terms_metadata3.txt",header = T, sep = "\t",row.names = 1))
rownames(term_sam)<-rownames(term_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")
term_sam$accession<-rownames(term_sam)

term_pseq

term_sam$outcome<-term_sam$outcome%>%
  str_replace_all("recovered", "Recovered")%>%
  str_replace_all("deceased","Deceased")%>%
  str_replace_all('stabilized',"Stabilized")
term_sam$sex<-term_sam$sex%>%
  str_replace_all("M", "male")%>%
  str_replace_all("F", "female")%>%
  str_replace_all("na", "<NA>") # this is mixing the string "<NA>" with actual NAs, probably not related to our problem, but def not a good idea in general

term_pseq <- phyloseq(term_counts_phy, term_tax_phy, sample_data(term_sam))
term_pseq# [ 14581 taxa and 167 samples ] [ 27077 taxa and 167 samples ]

filtme<-c("GO:0003674")
term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
filtme<-c("GO:0008150")
term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
term_pseq #[ 14579 taxa and 167 samples ]

term_pseq_no_neg<-term_pseq
term_pseq_no_neg<-subset_samples(term_pseq, sample_type!="neg_control")
term_pseq_no_neg# [ 14579 taxa and 162 samples ]:
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, sample_type!="Unknown")
term_pseq_no_neg#  [ 14579 taxa and 141 samples ]:
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, case!="Control_Sick")
term_pseq_no_neg# [ 14597 taxa and 105 samples ]
term_pseq_no_neg<-subset_samples(term_pseq_no_neg,publication!="Michalovich")
term_pseq_no_neg# [ 14597 taxa and 102 samples ]
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, bioproject!="PRJNA605907")
term_pseq_no_neg# [ 14597 taxa and 86 samples ]
term_pseq_no_neg<-prune_taxa(taxa = taxa_sums(term_pseq_no_neg)>0,x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ]
term_pseq_no_neg<-prune_samples(samples = sample_sums(term_pseq_no_neg)>0,x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ] # [ 25426 taxa and 86 samples ]
term_pseq_no_neg_gonames<-term_pseq_no_neg

meta<-meta(term_pseq_no_neg)%>%filter(case=="COVID19")
colnames(meta)

tally(sample_name~isolate_name,meta)

tax<-data.frame(tax_table(term_pseq_no_neg))
names<-paste(rownames(tax),tax$name,sep="-")
taxa_names(term_pseq_no_neg)<-names

save.image(file = "./images/0_preprocessing.RDA")

taxa_names(term_pseq_no_neg)

load(file = "./images/0_preprocessing.RDA")

library(mosaic)
library(microbiome)


meta<-as_tibble(meta(term_pseq))

tally(outcome~case=="COVID19",meta)

tally(case~publication=="Zhou",meta)

tally(publication~sample_name,meta)




