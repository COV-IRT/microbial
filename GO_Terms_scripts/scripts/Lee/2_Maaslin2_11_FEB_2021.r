library(tidyverse)
library(phyloseq)
library(Maaslin2)

setwd('/home/jovyan/work/microbial/GO_term_analysis/16_DEC_2020_dmm_troubleshooting_Lee/')

load(file = "././images/1_dmm_03222021.RDA")

term_pseq_no_neg_comp<-microbiome::transform(x = term_pseq_no_neg,transform = "compositional")

library(microbiome)

df_input_data<-data.frame(t(otu_table(term_pseq_no_neg_comp)))
df_input_metadata<-data.frame(sample_data(term_pseq_no_neg_comp))

write.table(x = df_input_metadata,file = "maaslin2_case_analysis_metadata.tsv",quote = F,sep = "\t",row.names = T,)

case_norm<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./results/terms_vs_case_comp_norm",
  min_abundance = 0.01, 
  min_prevalence = 0.1, 
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("sample_name","publication"),
  fixed_effects = c("case"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num,
 reference="case,COVID19")

term_pseq_outcome<-subset_samples(physeq = term_pseq_no_neg,outcome!="NA")
x<-sample_sums(term_pseq_outcome)
y<-as_tibble(x,,rownames="sample")%>%arrange(value)
meta<-as_tibble(meta(term_pseq_outcome),rownames="sample")

y<-full_join(y,meta)
library(ggpubr)
gghistogram(data = y,x = "value",y = "..count..",color = "publication",fill = "outcome",bins = 50)+xscale("log10")
# library(mosaic)
 favstats(data = y,x = y$value)

meta<-meta(term_pseq_outcome)
meta$outcome<-gsub("Stabilized","Survived",meta$outcome)
meta$outcome<-gsub("Recovered","Survived",meta$outcome)
sample_data(term_pseq_outcome)<-sample_data(meta)

sample_data(term_pseq_outcome)

remove_zhou<-c("SRR11092058","SRR11092057","SRR11092056","SRR11092064")

keep<-c("SRR11092061","CRR125941","CRR125940","CRR125939","CRR125938","CRR125937","CRR125936",
"CRR125935","SRR10903402","SRR10903401","SRR11092063","SRR11092062","SRR11092060","SRR11092059")


a<-as.data.frame(taxa_sums(term_pseq_outcome))

library(mosaic)

favstats(a$`taxa_sums(term_pseq_outcome`)

term_pseq_outcome_core<-prune_taxa(taxa_sums(term_pseq_outcome)>10,term_pseq_outcome)
term_pseq_outcome_core<-subset_samples(physeq = term_pseq_outcome_core,reads=="paired")
term_pseq_outcome_core<-prune_samples(x = term_pseq_outcome_core,samples = keep)
term_pseq_outcome
term_pseq_outcome_core

term_pseq_outcome_comp<-microbiome::transform(x = term_pseq_outcome_core,transform = "compositional")

COVID19_df_input_data<-data.frame(t(otu_table(term_pseq_outcome_comp)))
COVID19_df_input_metadata<-data.frame(sample_data(term_pseq_outcome_comp))

tally(~sample_name+publication,data=COVID19_df_input_metadata,"count")

tally(~sample_name+outcome,data=COVID19_df_input_metadata,"count")

tally(~outcome+publication,data=COVID19_df_input_metadata,"count")

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
#  random_effects = c("publication"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("publication"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c('sample_name'),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c('publication'),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)



subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c('reads'),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

COVID19_df_input_metadata_test<-data.frame(sample_data(term_pseq_outcome_comp))
COVID19_df_input_metadata_test$control<-as_factor(paste0(COVID19_df_input_metadata_test$publication,
    "_",
    COVID19_df_input_metadata_test$sample_name,
    "_",
    COVID19_df_input_metadata_test$reads))
COVID19_df_input_metadata_test$control

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("control"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome_res<-subset_outcome$results%>%filter(pval<=0.01)
subset_outcome_res

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("publication","sample_name"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome_res<-subset_outcome$results%>%filter(pval<=0.01)
dim(subset_outcome_res)

subset_outcome<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata_test,
  output="./results/subset_outcome",
  min_abundance = 0.00001, 
  min_prevalence = 0.01,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("publication","reads"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)

subset_outcome_res<-subset_outcome$results%>%filter(pval<=0.01)
dim(subset_outcome_res)

subset_outcome_res



tally(~sample_name+publication,data=COVID19_df_input_metadata_test,"count")

tally(~sample_name+outcome,data=COVID19_df_input_metadata_test)

tally(~sample_name+reads+publication,data=COVID19_df_input_metadata_test)

tally(~sample_name+reads+publication,data=COVID19_df_input_metadata_test)

tally(~control+outcome,data=COVID19_df_input_metadata_test)







tally(~sample_name+reads,meta,"count")

res<-case_norm$results%>%filter(pval<=0.05)

write.table(subset_outcome_res,"subset_outcome_res.tsv", sep="\t", row.names=F)

dim(subset_outcome_res)

Terms<-gsub("GO.","GO:",res$feature)
Terms<-gsub("[.]"," ",Terms)
Terms<-sub(" ","-",Terms)
Terms<-as_tibble(Terms)
Terms<-separate(data = Terms,col = value,sep = "-",into =  c("Term", "name"))

outcome_Terms<-gsub("GO.","GO:",subset_outcome_res$feature)
outcome_Terms<-gsub("[.]"," ",outcome_Terms)
outcome_Terms<-sub(" ","-",outcome_Terms)
outcome_Terms<-as_tibble(outcome_Terms)
outcome_Terms<-separate(data = outcome_Terms,col = value,sep = "-",into =  c("Term", "name"))

outcome_term_pseq_prune <- prune_taxa(taxa = outcome_Terms$Term,x =term_pseq_no_neg_gonames)
outcome_term_pseq_prune #[85 taxa and 141 samples ]
tax<-data.frame(tax_table(outcome_term_pseq_prune))
names<-paste(rownames(tax),tax$name,sep="-")
taxa_names(outcome_term_pseq_prune)<-names

term_pseq_prune <- prune_taxa(taxa = Terms$Term,x =term_pseq_no_neg_gonames)
term_pseq_prune #[85 taxa and 141 samples ]
tax<-data.frame(tax_table(term_pseq_prune))
names<-paste(rownames(tax),tax$name,sep="-")
taxa_names(term_pseq_prune)<-names

save.image("./images/2_maaslin2.rda")

write.table(tax,"subset_outcome_tax.tsv",sep = "\t",row.names = T,quote = F)

boxplot_abundance()

subset_outcome_res$feature


