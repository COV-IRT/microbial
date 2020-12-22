#if (!requireNamespace("BiocManager", quietly = TRUE))
#BiocManager::install(c("Maaslin2", "DESeq2",'phyloseq','microbiome','DirichletMultinomial','GenomicRanges'))

#install.packages("remotes")
#remotes::install_github("mikemc/speedyseq")
#install.packages(c('circlize','ggpubr','viridis','mosaic'))

library(tidyverse)
library(phyloseq)
library(microbiome)
library(DESeq2)
library(Maaslin2)
library(parallel)
library(DirichletMultinomial)
library(pheatmap)
library(ggpubr)
library(viridis)
library(mosaic)

#setwd('/home/jovyan/work/Aagaard_Raid3/microbial/GO_term_analysis/')
setwd('/media/jochum00/Aagaard_Raid3/jupyter_notebooks/jochum00_jupyter/microbial/GO_term_analysis/')
getwd()

ann_colors_old = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF","5"="#008280FF","6"="#BB0021FF","7"="#5F559BFF"),
  Publication=c("Huang"="dodgerblue4","Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Sample_Type =c("COVID_19"="firebrick",
                 "Healthy"="forestgreen",
                 "Community_acquired_pneumonia"="darkorange1",
                 "Obese_Asthma"="goldenrod3",
                 "Obese_Smoker"="goldenrod4", 
                 "Obese"="goldenrod1",
                 "Asthma"= "dodgerblue2", 
                 "Asthma_Smoker"="dodgerblue4",
                 "Asthma_Ex_smoker"="dodgerblue3",
                 "Smoker"="gray27",
                 "Obese_Asthma_Smoker"="black"),
  Case=c("Control_Neg"="grey75","Control_Unknown"="grey50","COVID19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Control_Healthy"="forestgreen","Control_Sick"="dodgerblue4"),
Outcome=c("Deceased"="black","Stabilized"="goldenrod4","Recovered"="forestgreen"))
ann_colors = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF"),#,"5"="#008280FF","6"="#BB0021FF","7"="#5F559BFF"),
  Publication=c("Michalovich"="#BB0021FF","Xiong"="#008280FF", "Shen"="#631879FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Case=c("COVID19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Control_Healthy"="forestgreen"),
 # Sample_Type=c("COVID_19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Healthy"="forestgreen"),
Outcome=c("Deceased"="black","Stabilized"="goldenrod4","Recovered"="forestgreen"))

raw<-as_tibble(read.table("../../Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))

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

bio_term<-bio%>%filter(type%in%c("bac","arc","vir"))%>%
    select(-type) %>%
    group_by(GO_term,namespace,depth,name)%>%
    summarise(across(.cols = where(is.numeric), sum))

mol_term<-mol%>%
    filter(type%in%c("bac","arc","vir"))%>%
    select(-type)%>%
    group_by(GO_term,namespace,depth,name)%>%
    summarise(across(.cols = where(is.numeric), sum))

bio_term_counts<-bio_term%>%
    select(-c(namespace,depth,name))
bio_term_tax<-bio_term%>%
    select(GO_term,namespace,depth,name)
mol_term_counts<-mol_term%>%
    select(-c(namespace,depth,name))
mol_term_tax<-mol_term%>%
    select(GO_term,namespace,depth,name)

bio_term_counts$namespace<-NULL
bio_term_counts$depth<-NULL
mol_term_counts$namespace<-NULL
mol_term_counts$depth<-NULL

# bio_term_tax<-bio_term%>%
#     select(GO_term,namespace,depth,name)
# mol_term_tax<-mol_term%>%
#     select(GO_term,namespace,depth,name)

bio_term_counts_df<-data.frame(bio_term_counts, row.names=1)
bio_term_tax_df<-data.frame(bio_term_tax, row.names=1)
mol_term_counts_df<-data.frame(mol_term_counts, row.names=1)
mol_term_tax_df<-data.frame(mol_term_tax, row.names=1)

dim(bio_term_counts)
dim(bio_term_tax)
dim(mol_term_counts)
dim(mol_term_tax)
#mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax), errorIfNULL = T)

bio_term_counts_phy <- otu_table(bio_term_counts_df, taxa_are_rows=TRUE)
bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax_df), errorIfNULL=TRUE)
mol_term_counts_phy<-otu_table(mol_term_counts_df, taxa_are_rows = T)
mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax_df), errorIfNULL = T)

getwd()

bio_term_sam<-as.data.frame(read.table("../../Combined_BALF_GO_Terms_metadata2.txt",header = T, sep = "\t",row.names = 1))

rownames(bio_term_sam)<-rownames(bio_term_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")
bio_term_sam$accession<-rownames(bio_term_sam)

bio_term_sam$outcome<-bio_term_sam$outcome%>%
str_replace_all("recovered", "Recovered")%>%
str_replace_all("deceased","Deceased")%>%
str_replace_all('stabilized',"Stabilized")
###DONT FORGET TO DELETE THESE LINES LATER AFTER YOUR DONE PLAYING AROUND
#str_replace_all('Stabilized',"Survived")%>%
#str_replace_all("Recovered", "Survived")
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

meta<-meta(term_pseq)
meta$idsums<-sample_sums(term_pseq)
meta<-as_tibble(meta)%>%arrange(desc(idsums))
gghistogram(data = meta,
            title = "GO Term Totals by Publication",
            legend="bottom",
            x="idsums",
            y="..count..",
            color="case",
            fill="case",
            facet="publication", 
            bins = 30,
            rug=T,
            add_density=T,
            alpha=0.3,
            add = c("mean"))+
scale_x_log10()+
scale_fill_discrete(type = ann_colors_old)+
scale_color_discrete(type = ann_colors_old)
tally(publication~case,meta,"count")

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

meta<-meta(term_pseq_no_neg)
tally(~case, data = meta, format ="count")
tally(~outcome, data = meta,format ="count")
tally(case~publication,data = meta,format ="count")

#names<-paste(taxa_names(term_pseq_no_neg),get_taxa_unique(term_pseq_no_neg,taxonomic.rank = "name" ),sep = "-")
#taxa_names(term_pseq_no_neg)<-names

tax<-data.frame(tax_table(term_pseq_no_neg))
names<-paste(rownames(tax),tax$name,sep="-")
length(names)
taxa_names(term_pseq_no_neg)<-names

meta<-meta(term_pseq_no_neg)
meta$idsums<-sample_sums(term_pseq_no_neg)
meta<-as_tibble(meta)%>%arrange(desc(idsums))
gghistogram(data = meta,
            title = "GO Term Totals by Case",
            legend="bottom",
            x="idsums",
            y="..count..",
            color="case",
            fill="case", 
            bins = 30,
            rug=T,
            add_density=T,
            alpha=0.3,
            add = c("mean"))+
scale_x_log10()+
scale_fill_discrete(type = ann_colors)+
scale_color_discrete(type = ann_colors)

meta<-meta(term_pseq_no_neg)
meta$idsums<-sample_sums(term_pseq_no_neg)
meta<-as_tibble(meta)%>%arrange(desc(idsums))
gghistogram(data = meta,
            title = "GO Term Totals by Publication",
            legend="bottom",
            x="idsums",
            y="..count..",
            color="publication",
            fill="publication", 
            bins = 30,
            rug=T,
            add_density=T,
            alpha=0.3,
            add = c("mean"))+
scale_x_log10()+
scale_fill_discrete(type = ann_colors)+
scale_color_discrete(type = ann_colors)

term_pseq_no_neg_comp<-microbiome::transform(x = term_pseq_no_neg,transform = "compositional")

sample_info_tab<-sample_data(term_pseq_no_neg)
sample_info_tab_phy <- sample_data(sample_info_tab)
deseq_counts<-phyloseq_to_deseq2(physeq = term_pseq_no_neg,design = ~ 1) 
deseq_counts_vst<-varianceStabilizingTransformation(deseq_counts)
#deseq_counts_vst <- estimateSizeFactors(deseq_counts, type = "poscounts")
vst_trans_count_tab <- assay(deseq_counts_vst)

vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_tax_phy <- tax_table(term_pseq_no_neg)
vst_physeq <- phyloseq(vst_count_phy, vst_tax_phy,sample_data(term_pseq_no_neg))
#vst_physeq_comp<-microbiome::transform(x = vst_physeq,transform = "compositional")

#dir.create("R_Maaslin2") # Create a new directory
setwd('/media/jochum00/Aagaard_Raid3/jupyter_notebooks/jochum00_jupyter/microbial/GO_term_analysis/R_Maaslin2/ag1')
#setwd("/home/jovyan/work/Jochum_3/jupyter_lab/GO_term_analysis/R_Maaslin2/ag1") # Change the current working directory 
getwd() #check if directory has been successfully changed

df_input_data<-data.frame(t(otu_table(term_pseq_no_neg_comp)))
df_input_metadata<-data.frame(sample_data(term_pseq_no_neg_comp))

df_input_data2<-data.frame(t(otu_table(vst_physeq)))
df_input_metadata2<-data.frame(sample_data(vst_physeq))

term_pseq_no_neg_comp#  [ 13534 taxa and 86 samples ]

num<-as.numeric(100)
case_norm<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./terms_vs_case_comp_norm",
  min_abundance = 0.01, #UPDATE for more terms: lowered miminum rel. abundance 0.3%
  min_prevalence = 0.1, 
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  #random_effects = c("sample_name","publication","collection_location","sequence_type"),
  random_effects = c("sample_name","publication"),
  fixed_effects = c("case"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num,
  reference=c("case,COVID19"))

#I need to figure out how to change how to pivot wider the case colum for an age analysis
df_input_metadata$age<-as.numeric(df_input_metadata$age)

age<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./terms_vs_case_comp_norm_age_fixed",
  min_abundance = 0.01,
  min_prevalence = 0.1,
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.25,
  #random_effects = c("sample_name","publication","collection_location","sequence_type"),
  random_effects = c("sample_name","publication"),
  fixed_effects = c("age","case"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num,
  reference=c("case,COVID19"))

term_pseq_outcome<-subset_samples(physeq = term_pseq_no_neg,outcome!="NA")
#term_pseq_outcome<-subset_samples(physeq = term_pseq_no_neg, case=="COVID19")
term_pseq_outcome# [ 13534 taxa and 25 samples ]
sample_info_tab<-sample_data(term_pseq_outcome)
sample_info_tab_phy <- sample_data(sample_info_tab)
deseq_counts<-phyloseq_to_deseq2(physeq = term_pseq_outcome,design = ~ 1)
deseq_counts_vst<-varianceStabilizingTransformation(deseq_counts)
#deseq_counts_vst <- estimateSizeFactors(deseq_counts, type = "poscounts")
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_tax_phy <- tax_table(term_pseq_no_neg)
vst_physeq <- phyloseq(vst_count_phy, vst_tax_phy,sample_data(term_pseq_outcome))
vst_physeq_comp<-microbiome::transform(x = vst_physeq,transform = "compositional")
#df_input_data<-data.frame(t(otu_table(vst_physeq_comp)))
#df_input_metadata<-data.frame(sample_data(vst_physeq_comp))
term_pseq_outcome_comp<-microbiome::transform(x = term_pseq_outcome,transform = "compositional")
df_input_data<-data.frame(t(otu_table(term_pseq_outcome_comp)))
df_input_metadata<-data.frame(sample_data(term_pseq_outcome_comp))

df_input_metadata$age<-as.numeric(df_input_metadata$age)

subset_outcome<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./subset_outcome_age_fixed",
  min_abundance = 0.001,
  min_prevalence = 0.05,
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.5,
  random_effects = c("sample_name","publication"),
  fixed_effects = c("age","outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = num,
  reference="outcome,Deceased")

#df_input_metadata$age<-as.numeric(df_input_metadata$age)
subset_outcome<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./subset_outcome",
  min_abundance = 0.001,
  min_prevalence = 0.05,
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.5,
  random_effects = c("sample_name","publication"),
  fixed_effects = c("outcome"),
  correction="BH",
  standardize = TRUE,
  cores = 48,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = num,
  reference="outcome,Deceased")

#case$results

res<-case_norm$results%>%filter(pval<=0.05)
res

Terms<-gsub("GO.","GO:",res$feature)
Terms<-gsub("[.]"," ",Terms)
Terms<-sub(" ","-",Terms)
Terms<-as_tibble(Terms)
Terms<-separate(data = Terms,col = value,sep = "-",into =  c("Term", "name"))

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

term_pseq_prune <- prune_taxa(taxa = Terms$Term,x =term_pseq_no_neg)
term_pseq_prune 
term_pseq_prune<-prune_taxa(taxa_sums(term_pseq_prune)>0,term_pseq_prune)
term_pseq_prune
term_pseq_prune<-prune_samples(sample_sums(term_pseq_prune)>0,term_pseq_prune)
term_pseq_prune<-prune_taxa(taxa_sums(term_pseq_prune)>0,term_pseq_prune)
term_pseq_prune# 35 taxa and 86 samples 

tax<-data.frame(tax_table(term_pseq_prune))
names<-paste(rownames(tax),tax$name,sep="-")
taxa_names(term_pseq_prune)<-names

dat <- abundances(term_pseq_prune)
count <- as.matrix(t(dat))

library(parallel)
fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)

#fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)

lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

#identify the number of clusters that best fits the model

best <- fit[[which.min(lplc)]]
#best <-fit[[4]]
best

#save.image(file = "go_terms_dmm.rdata")

options(repr.plot.width=15, repr.plot.height=15)
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =log2, lblwidth = 0.2 * nrow(count))

options(repr.plot.width=15, repr.plot.height=15)
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =sqrt, lblwidth = 0.2 * nrow(count))

mixturewt(best)

write.table(fitted(best),"GO_TERMS_DMM_contributions.tsv", sep="\t")

ass <- apply(mixture(best), 1, which.max)
write.table(ass,"GO_TERMS_DMM_groups.tsv",sep="")

#add the dmm group to the metadata
sample_data(term_pseq_prune)$dmn<-ass
term_pseq_prune_comp<-microbiome::transform(term_pseq_prune,"compositional")
#melt the phyloseq object into tidy form
tmp<-psmelt(term_pseq_prune_comp)
tmp<-as_tibble(tmp)

options(repr.plot.width=5, repr.plot.height=5)
gghistogram(tmp,x = "Abundance",y = "..count..")+scale_x_log10() #move each dmm group into a colum of its own

#tmp$log2Abundance<-log2(tmp$Abundance)

#subset the dataset to only include the case, Go_term, count, and dmm group.
#obtain the avergage count for each Go term
#order the go terms from hight to lowest count

d2<-tmp %>%
  select(case,OTU,Abundance, dmn)%>%
  group_by(OTU,case, dmn) %>%
  summarise(avg = mean(Abundance)) %>%
  arrange(desc(avg))

#d2$avg<-sqrt(d2$avg)

#get the total count of the go terms and oder from greates to lowest


d3<-tidyr::spread(d2,dmn,avg)
d3$tot<-rowSums(d3[3:6], na.rm = T)
d3<-d3%>%arrange(desc(tot))
d3$tot<-NULL
head(d3)

d3<-d3%>%gather(data = d3,avg,3:6)
colnames(d3)<-c("name","case", "dmn","avg")
d3

d4<-d3%>%filter(avg>0.01)%>%arrange(name,case,dmn)
#d4<-d4[1:108,]
dim(d4)
dim(d3)

my_pal<-viridis(n = 256, alpha = 1, begin = 0, end = 1, direction = 1)

options(repr.plot.width=10, repr.plot.height=8)
a<-ggplot(data = d3,mapping = aes(x = factor(dmn),y =reorder(name,avg),size=avg,color=avg))+
geom_point()+
#theme(text=element_text(size=20))+
scale_colour_gradientn(colours = my_pal)+
facet_grid(facets = ~ case)+
theme_minimal()

a

options(repr.plot.width=10, repr.plot.height=8)  
ggballoonplot(d3, y ="name",x = "dmn", size = "avg", facet.by = "case",fill = "avg",ggtheme = theme_minimal())+
     guides(size = FALSE)+
    font("y.text", size = 12)+scale_fill_viridis_c()

save.image("GO_TERM_Maaslin2_29_NOV_2020.rda")
#load.Rdata("GO_TERM_Maaslin2_19_NOV_2020.rda")

sam<-as_tibble(sample_data(term_pseq_prune))
tbl<-tally(x = case~dmn,sam, format = "count")

case<-as.vector(tally(~case,sam))
case
dmn<-as.vector(tally(~dmn,sam))
dmn
cor.test(case,dmn)
cor(case,dmn)
cov(case,dmn)

out_tbl<-tally(x = case~dmn,sam)
result2 <-chisq.test(table(out_tbl))
out_tbl
result2
chisq(result2)



###########################################
###Dont forget to save you shit HERE#######
###########################################
#save.image(file = "term_go_terms_dmm.rdata")
#

count<-abundances(term_pseq_prune)
select <- order(rowMeans(count),decreasing=TRUE)
select3<-order(rowSdDiffs(count),decreasing=T)
select2<-log1p((count)[select3,])
tmp<-rownames(select2)
dim(select2)

#select2$mol<-tmp
select_tibb<-as_tibble(select2)
rownames(select_tibb)<-rownames(select2)
write.table(select2,file = "GO_Terms_results")

#library(matrixTests)
#library(genefilter)
#select3<-as_tibble(select2)%>%summarise(std=rowFtests(select2))%>% arrange(desc(std))
#select2$mol<-tmp
select2<-as_tibble(select2)
#select2$mean<-rowMeans(select2)

rownames(select2)<-tmp
sam<-data.frame(sample_data(term_pseq_prune))
df<-as.data.frame(sample_data(term_pseq_prune))
df<-as_tibble(df)
#df<-df%>%select(dmn, publication, sample_type,case,outcome)#dmn,body_site)
#update, removing sample type column after removing sick and micalovich samples
df<-df%>%select(dmn, publication, case,outcome)#dmn,body_site)
df<-as.data.frame(df)

#select2<-sqrt((count)[select,])
dim(select2)
colnames(select2) <- colnames(otu_table(term_pseq_prune))
length(rownames(select2))
length(row.names(df))
row.names(df) <- colnames(select2)

library(ggsci)

mypal <- pal_aaas("default", alpha = 1)(10)
mypal
library("scales")
library(RColorBrewer)
library(viridis)
library(pheatmap)
df_row<-as.data.frame(fitted(best))
#colnames(df_row)<-c("1","2","3","4","5","6")
df<-data.frame(df)
#update, removing sampleType column after removinch "sick" samples
#colnames(df)<-c("dmm_cluster", "Publication","Sample_Type","Case","Outcome")
colnames(df)<-c("dmm_cluster", "Publication","Case","Outcome")

# Specify colors

ann_colors = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF","5"="#008280FF","6"="#BB0021FF","7"="#5F559BFF"),
  Publication=c("Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Sample_Type =c("COVID_19"="firebrick",
                 "Healthy"="forestgreen",
                 "Community_acquired_pneumonia"="darkorange1",
                 "Obese_Asthma"="goldenrod3",
                 "Obese_Smoker"="goldenrod4", 
                 "Obese"="goldenrod1",
                 "Asthma"= "dodgerblue2", 
                 "Asthma_Smoker"="dodgerblue4",
                 "Asthma_Ex_smoker"="dodgerblue3",
                 "Smoker"="gray27",
                 "Obese_Asthma_Smoker"="black"),
  Case=c("COVID19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Control_Healthy"="forestgreen","Control_Sick"="dodgerblue4"),
Outcome=c("Deceased"="black","Stabilized"="goldenrod4","Recovered"="forestgreen"))

ann_colors = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF"),#,"5"="#008280FF","6"="#BB0021FF","7"="#5F559BFF"),
  Publication=c("Xiong"="#008280FF", "Shen"="#631879FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Case=c("COVID19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Control_Healthy"="forestgreen"),
 # Sample_Type=c("COVID_19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Healthy"="forestgreen"),
Outcome=c("Deceased"="black","Stabilized"="goldenrod4","Recovered"="forestgreen"))

select3<-select2%>%filter(rowMeans(select2)>8)
dim(select2)
dim(select3)

options(repr.plot.width=20, repr.plot.height=8)
xx <- pheatmap(mat = select3,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),
               annotation_col=df,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
              annotation_row = df_row)

options(repr.plot.width=24, repr.plot.height=5)
xx <- pheatmap(mat = select2,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10),
               annotation_col=df,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               annotation_row = df_row,
               cluster_row = T,
               cluster_col=T)

#BiocManager::install(c('FactoMineR','factoextra'))

library('FactoMineR')
library('factoextra')
options(repr.plot.width=5, repr.plot.height=5)
a<-tally(~case+dmn,meta(term_pseq_prune))
res.CA<-CA(a,graph=F)
#res.CA
fviz_ca_biplot(res.CA,
               repel=T,
               col.col="cos2",
               col.row="cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

library('FactoMineR')
library('factoextra')
options(repr.plot.width=8, repr.plot.height=8)
a<-tally(~outcome+dmn,meta(term_pseq_prune))
res.CA<-CA(a,graph=F)
#res.CA
fviz_ca_biplot(res.CA,
               repel=T,
               col.col="cos2",
               col.row="cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

##################ANYTHING BELOW THIS LINE IS TRASH FOR NOW###########
####################################################
# Correlation between dmm groups and cases
#################################################




library(matrixTests)
p<-psmelt(term_pseq_prune)
dmn_sum<-p%>%
  select(Sample,OTU,dmn,Abundance)%>%
  group_by(OTU,dmn)%>%
  pivot_wider(id_cols = c(Sample,dmn), names_from = OTU,values_from = Abundance)
write.table(dmn_sum, "dmn_sum.tsv",sep="\t")
# an example with offsets from Venables & Ripley (2002, p.189)
dmn_case<-pairwise.wilcox.test(p$Abundance, p$case,p.adjust.method = "BH")
dmn_wilcox<-pairwise.wilcox.test(p$Abundance, p$dmn,p.adjust.method = "BH")


p$lg<-log1p(p$Abundance)
p$dmn<-as_factor(p$dmn)
p$case<-factor(x = p$case, levels = c("Control_Healthy","Control_Sick","COVID19"))


library(lmerTest)
library(lmer)
case_dmn_glm <- lmer(lg ~ dmn +(dmn|case),data = p)
f<-ls_means(case_dmn_glm)

# b<-glm( lg ~ dmn+case, data = p, family = gaussian)
library(mosaic)
msummary(b)
b$coefficients
c<-anova(b)
d<-aov(b)
e<-TukeyHSD(d)
library(lmerTest)
f<-ls_means(case_dmn_glm)

msummary(e$case)
e$case
mat<-dmn_sum[,3:length(dmn_sum)]
krus_dmn<-col_kruskalwallis(x =mat,g = dmn_sum$dmn)
welch_dmn<-col_oneway_welch(x=mat, g=dmn_sum$dmn)

krus_dmn_sig<-krus_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_dmn_sig<-welch_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)

TukeyHSD(welch_dmn_sig)
library(mosaic)
msummary(welch_dmn_sig)
dim(krus_dmn_sig)
dim(welch_dmn_sig)
krus_sig
write.table(krus,"krus_dmn.tsv",sep = "\t")

colnames(krus_sig)
colnames(krus_sig)
krus_sig_names<-intersect(rownames(krus_case_sig), rownames(krus_dmn_sig))
welch_sig_names<-intersect(rownames(welch_case_sig), rownames(welch_dmn_sig))
sig_names<-intersect(krus_sig_names,welch_sig_names)
