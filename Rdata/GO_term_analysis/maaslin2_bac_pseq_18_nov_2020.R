# library(phyloseq)
# library(microbiome)
# library(DirichletMultinomial)
# library(reshape2)
# library(magrittr)
# library(dplyr)
# library(tidyr)
# library(pheatmap)
# library(ggpubr)
# library(ggsci)
######################################################
#    PART 1: DATA WRANGLING
######################################################
"###################MIKE LEE PLEASE HELPPPP#############################
########WHY ARE THE GO_TERM NAMES BREAKING THE DATAFRAMES!?!?!#########"

#load libraries
library(phyloseq)
library(tidyverse)


raw<-as_tibble(read.table("Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
# A tibble: 47,233 x 2,020     # good so far now

#do a little regex and fix some stuff
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

#transform the raw table by type of count (euk, term, bac, arc)
df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","Bacteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)


####There are multiple processes and values for a single sample so you cant convert the sample to columns
#make individual tibbles for biological processes and molecular fxn
bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")

#make individual tibbles for each type (bac, euk, term, arc, vir, etc)
bio_bac<-bio%>%filter(type=="bac")%>%select(-type)
bio_term<-bio%>%filter(type=="term")%>%select(-type)
mol_bac<-mol%>%filter(type=="bac")%>%select(-type)
mol_term<-mol%>%filter(type=="term")%>%select(-type)

#subselect tibbles for only the counts and go terminology
bio_bac_counts<-bio_bac%>%select(-c(namespace,depth,name))
bio_bac_tax<-bio_bac%>%select(GO_term,namespace,depth,name)
mol_bac_counts<-mol_bac%>%select(-c(namespace,depth,name))
mol_bac_tax<-mol_bac%>%select(GO_term,namespace,depth,name)

#convert them to dataframes for downstream import to phylsoeq
bio_bac_counts<-data.frame(bio_bac_counts, row.names=1)
bio_bac_tax<-data.frame(bio_bac_tax, row.names=1)
mol_bac_counts<-data.frame(mol_bac_counts, row.names=1)
mol_bac_tax<-data.frame(mol_bac_tax, row.names=1)

#convert the dataframes into phyloseq formats
bio_bac_counts_phy <- otu_table(bio_bac_counts, taxa_are_rows=TRUE)
bio_bac_tax_phy <- tax_table(as.matrix(bio_bac_tax), errorIfNULL=TRUE)
mol_bac_counts_phy<-otu_table(mol_bac_counts, taxa_are_rows = T)
mol_bac_tax_phy<-tax_table(as.matrix(mol_bac_tax), errorIfNULL = T)


#import your metadata
bio_bac_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata.txt",header = T, sep = "\t",row.names = 1))
#a little regex to fix the stupid filename
rownames(bio_bac_sam)<-rownames(bio_bac_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")
bio_bac_sam$accession<-rownames(bio_bac_sam)

# making physeq object
bio_bac_pseq <- phyloseq(bio_bac_counts_phy, bio_bac_tax_phy, sample_data(bio_bac_sam))
mol_bac_pseq<-phyloseq(mol_bac_counts_phy,mol_bac_tax_phy, sample_data(bio_bac_sam))
library(speedyseq)
bac_pseq<-merge_phyloseq(bio_bac_pseq,mol_bac_pseq)
bac_pseq #[ 13846 taxa and 167 samples ]:
#filter out the negative control and unknown samples
bac_pseq_no_neg<-subset_samples(bac_pseq, sample_type!="neg_control")
bac_pseq_no_neg# [ 13846 taxa and 162 samples ]:
bac_pseq_no_neg<-subset_samples(bac_pseq_no_neg, sample_type!="Unknown")
bac_pseq_no_neg# [ 13846 taxa and 141 samples ]:
