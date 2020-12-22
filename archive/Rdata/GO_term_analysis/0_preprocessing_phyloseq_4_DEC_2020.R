library(tidyverse)
library(phyloseq)
setwd("d:/github/microbial/Rdata/GO_term_analysis/")
getwd()

raw<-as_tibble(read.table("./Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","termteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)

df$depth<-as.character(df$depth)
#SIDE NOTE:There are multiple processes and values for a single sample so you cant convert the sample to columns

bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")


#REMEMBER TO FIX THIS AGAIN LATER 4 DEC 2020
# TROUBLESHOOTING CODE ONLY
#######################################################

####There are multiple processes and values for a single sample so you cant convert the sample to columns
  #Old WAY
    
  #make individual tibbles for each type (bac, euk, term, arc, vir, etc)
  bio_term<-bio%>%filter(type=="bac")%>%select(-type)
  mol_term<-mol%>%filter(type=="bac")%>%select(-type)
    
  #subselect tibbles for only the counts and go terminology
  bio_term_counts<-bio_term%>%select(-c(namespace,depth,name))
  bio_term_tax<-bio_term%>%select(GO_term,namespace,depth,name)
  mol_term_counts<-mol_term%>%select(-c(namespace,depth,name))
  mol_term_tax<-mol_term%>%select(GO_term,namespace,depth,name)
  
  #convert them to dataframes for downstream import to phylsoeq
  bio_term_counts<-data.frame(bio_term_counts, row.names=1)
  bio_term_tax<-data.frame(bio_term_tax, row.names=1)
  mol_term_counts<-data.frame(mol_term_counts, row.names=1)
  mol_term_tax<-data.frame(mol_term_tax, row.names=1)
  
  #convert the dataframes into phyloseq formats
  bio_term_counts_phy <- otu_table(bio_term_counts, taxa_are_rows=TRUE)
  bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax), errorIfNULL=TRUE)
  mol_term_counts_phy<-otu_table(mol_term_counts, taxa_are_rows = T)
  mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax), errorIfNULL = T)

####################################################
# COMMENTED OUT ORIGINAL CODE
#############################################
# bio_term<-bio%>%filter(type%in%c("bac","arc","vir"))%>%
#   select(-type) %>%
#   group_by(GO_term,namespace,depth,name)%>%
#   summarise(across(.cols = where(is.numeric), sum))
# 
# mol_term<-mol%>%
#   filter(type%in%c("bac","arc","vir"))%>%
#   select(-type)%>%
#   group_by(GO_term,namespace,depth,name)%>%
#   summarise(across(.cols = where(is.numeric), sum)) 
# bio_term_counts<-bio_term%>%
#   select(-c(namespace,depth,name))
# bio_term_tax<-bio_term%>%
#   select(GO_term,namespace,depth,name)
# mol_term_counts<-mol_term%>%
#   select(-c(namespace,depth,name))
# mol_term_tax<-mol_term%>%
#   select(GO_term,namespace,depth,name)

# bio_term_counts$namespace<-NULL
# bio_term_counts$depth<-NULL
# mol_term_counts$namespace<-NULL
# mol_term_counts$depth<-NULL

# bio_term_counts_df<-data.frame(bio_term_counts, row.names=1)
# bio_term_tax_df<-data.frame(bio_term_tax, row.names=1)
# mol_term_counts_df<-data.frame(mol_term_counts, row.names=1)
# mol_term_tax_df<-data.frame(mol_term_tax, row.names=1)
# 
# dim(bio_term_counts)
# dim(bio_term_tax)
# dim(mol_term_counts)
# dim(mol_term_tax)
# #mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax), errorIfNULL = T)
# 
# bio_term_counts_phy <- otu_table(bio_term_counts_df, taxa_are_rows=TRUE)
# bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax_df), errorIfNULL=TRUE)
# mol_term_counts_phy<-otu_table(mol_term_counts_df, taxa_are_rows = T)
# mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax_df), errorIfNULL = T)
#################################################################


bio_term_sam<-as.data.frame(read.table("./Combined_BALF_GO_Terms_metadata2.txt",header = T, sep = "\t",row.names = 1))

rownames(bio_term_sam)<-rownames(bio_term_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")
bio_term_sam$accession<-rownames(bio_term_sam)

bio_term_sam$outcome<-bio_term_sam$outcome%>%
  str_replace_all("recovered", "Recovered")%>%
  str_replace_all("deceased","Deceased")%>%
  str_replace_all('stabilized',"Stabilized")

###############################################

bio_term_sam$sex<-bio_term_sam$sex%>%
  str_replace_all("M", "male")%>%
  str_replace_all("F", "female")%>%
  str_replace_all("na", "<NA>")

bio_term_pseq <- phyloseq(bio_term_counts_phy, bio_term_tax_phy, sample_data(bio_term_sam))
mol_term_pseq<-phyloseq(mol_term_counts_phy,mol_term_tax_phy, sample_data(bio_term_sam))
term_pseq<-merge_phyloseq(bio_term_pseq,mol_term_pseq)
term_pseq# [ 14581 taxa and 167 samples ]
filtme<-c("GO:0003674")
term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
filtme<-c("GO:0008150")
term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
term_pseq

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
term_pseq_no_neg# [ 13534 taxa and 86 samples ]