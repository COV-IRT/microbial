library(phyloseq)
library(tidyverse)
#setwd("/media/jochum00/Aagaard_Raid3/microbial/step2_kraken2_analysis/")
setwd("c:/github/microbial/step2_kraken2_analysis/")

#filter out only the 86 samples used in the GO term analysis
sam<-as_tibble(read.table("Combined_BALF_GO_Terms_metadata.tsv",header = T,sep = "\t"))
go_term_samples<-as_tibble(read.table("GO_TERM_samples.tsv",header=T,sep="\t"))
negs<-as_tibble(c("CRR125995",'CRR125996',"CRR125997","CRR125998"))
negs$accession<-negs$value
negs$value<-NULL
go_term_samples
keep<-bind_rows(go_term_samples,negs)
keep
sam<-left_join(keep,sam)
sam
#make a vector of the sample names
vec<-sam$accession
#filter out only the 86 samples used in the go term analysis
counts<-as_tibble(read.table("counts.tsv",header=T,sep="\t"))%>%
  pivot_longer(cols = -c(taxid),names_to="accession")%>%
  filter(accession %in% vec)%>%
  pivot_wider(names_from=accession)
#import the taxonomy and merge it up with the counts
tax<-as_tibble(read.table("taxonomy.tsv",header=T,sep="\t"))#,rownames = "taxid")
# convert to df
counts_df<-data.frame(counts,row.names = 1)
tax_df<-data.frame(tax,row.names = 1)
sam_df<-data.frame(sam,row.names=1)
#make a phyloseq obj
count_pseq<-otu_table(count_df,taxa_are_rows = T)
tax_pseq<-tax_table(as.matrix(tax_df),errorIfNULL = T)
sam_pseq<-sample_data(sam_df)
pseq<-phyloseq(count_pseq,tax_pseq,sam_pseq)
pseq<-subset_taxa(physeq = pseq,domain=="Bacteria")
pseq
# Detection threshold 0 (strictly greater by default);
# Prevalence threshold 50 percent (strictly greater by default)
#pseq <- core(pseq, 0, 50/100)
pseq



########################################################################################
library(decontam)
#export the shen et al. samples for parsing with decontam
pseq_shen<-subset_samples(pseq,bioproject=="CRA002476")
pseq_not_shen<-subset_samples(pseq,bioproject!="CRA002476")

count_df_shen<-abundances(pseq_shen)
tax_table(pseq_shen) %>% as("matrix")# %>% subset(select = colSums(!is.na(.)) > 0) %>% as_tibble(rownames = "OTU") -> taxon.tbl

tax_df_shen<-as.matrix(tax_table(pseq_shen))
tax_df
tax_df_shen<-data.frame(tax_table(pseq_shen))
tax_df_shen
sam_shen<-meta(pseq_shen)
vector_for_decontam <- sam_shen$case=="Control_Neg"
vector_for_decontam                    


count_df_shen<-abundances(pseq_shen)
colSums(count_df_shen)

contam_df_shen <- isContaminant(t(count_df_shen), neg=vector_for_decontam)


table(contam_df_shen$contaminant) # identified 4644 as contaminants 




# getting vector holding the identified contaminant IDs
contam_shen <- row.names(contam_df_shen[contam_df_shen$contaminant == TRUE, ])

contam_shen<-as.factor(contam_shen)
contam_shen



count_df_shen<-as_tibble(count_df_shen,rownames = "taxid")

count_df_shen[1]<-NULL

clean_count_df_shen<-count_df_shen%>%filter(!taxid%in%contam_shen)

a<-clean_count_df_shen%>%select(-taxid)
b<-count_df_shen%>%select(-taxid)
perc_remain<-(colSums(a)/colSums(b))*100



sam_shen
contaminants<-tax%>%filter(taxid %in% contam_shen)


not_contaminants<-tax%>%filter(!taxid %in% contam_shen)

contaminants
not_contaminants



not_contaminants
write.table(x = contaminants,file = "Shen_contaminants.tsv",sep = "\t",row.names = F,quote = F)
write.table(x = not_contaminants,file = "Shen_not_contaminants.tsv",sep = "\t",row.names = F,quote = F)

contam_shen<-as.factor(contam_shen)

contaminants$taxid<-as.character(contaminants$taxid)
not_contaminants$taxid<-as.character(not_contaminants$taxid)

pseq_shen_decontam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)
pseq_shen_decontam


library(microbiome)
meta(pseq_shen_contam)$case

pseq_shen_contam<-prune_taxa(taxa = not_contaminants$taxid,pseq_shen)

pseq_shen_decontam
decontam_melt<-psmelt(pseq_shen_decontam)
contam_melt<-psmelt(pseq_shen_contam)

decontam_melt<-as_tibble(decontam_melt,rownames = "taxid")
contam_melt<-as_tibble(contam_melt,rownames = "taxid")

decontam_melt_otu<-decontam_melt%>%select(Abundance,taxid,domain,phylum,class,order,family,genus,species,OTU,sample_name,sample_type)
decontam_melt_otu$contam<-FALSE
contam_melt_otu<-contam_melt%>%select(Abundance,taxid,domain,phylum,class,order,family,genus,species,OTU,sample_name,sample_type)
contam_melt_otu$contam<-TRUE


melt<-full_join(decontam_melt_otu,contam_melt_otu)
melt




d<-melt%>%
  select(Abundance,taxid,domain,phylum,class,order,family,genus,species,sample_type,sample_name,contam)%>%arrange(sample_type)

d

tally(taxid+contam~count(Abundance),
      data = d, margins = TRUE)









d$
library(ggpubr)
ggboxplot(d,x = "sample_type",y = "Abundance",palette = "aaas",facet.by = "phylum")+
  yscale("log10",.format = T)+rotate_x_text()

n_distinct(d$family)
library(mosaic)
range(d$Abundance)
d$sample_name
d%>%tally(~contam+Abundance)
d$sample_type
tally(sample_type~count(Abundance),
      data = d, margins = TRUE)
library(mosaic)
ggbarplot(data = d,x = "sample_type",y = "Abundance",color = "contam",palette = "npg",)


favstats(x =Abundance~contam+sample_name,data = d)
gghistogram(data = d,x = "Abundance",y = "..density..",color = "sample_type",palette = "aaas",facet.by = "sample_type")
