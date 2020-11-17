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
#remove duplicates samples by sample name
bac_pseq_no_neg_unique<-prune_samples(!duplicated(sample_data(bac_pseq_no_neg)$sample_name), x =  bac_pseq_no_neg)
#bac_pseq_no_neg_unique #[ 13846 taxa and 126 samples ]:
bac_pseq_no_neg_unique<-subset_samples(bac_pseq_no_neg_unique, publication!="Michalovich")
bac_pseq_no_neg_unique

#filter out empty GO Terms and empty samples
bac_pseq_prune <- prune_samples(sample_sums(bac_pseq_no_neg_unique) > 1, bac_pseq_no_neg_unique)
bac_pseq_prune# [ 13846 taxa and 126 samples ]:
bac_pseq_prune <- prune_taxa(taxa_sums(bac_pseq_prune) > 1, bac_pseq_prune)
bac_pseq_prune #[ 13476 taxa and 126 samples ]:
##########################################################
#Identify case associated significant taxa
#####################################################
library(DESeq2)
sample_info_tab<-sample_data(bac_pseq_prune)
sample_info_tab_phy <- sample_data(sample_info_tab)
deseq_counts<-phyloseq_to_deseq2(physeq = bac_pseq_prune,design = ~ 1) 
deseq_counts_vst <- estimateSizeFactors(deseq_counts, type = "poscounts")
vst_trans_count_tab <- assay(deseq_counts_vst)

#YAAAAAAAAAAAAAAAASSSSSSSSSSSSS
vst_trans_count_tab2 <- limma::removeBatchEffect(vst_trans_count_tab, sample_info_tab$publication)
#IT FIXED THE BATCH EFFECT WORRRRRRRRRRRRRRRRRRRKED##############

####lets see if there are still batch effect stuff going on
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
plot(euc_clust) 
euc_dend <- as.dendrogram(euc_clust, hang=-15)
library(dendextend)
dend_cols <- sample_info_tab$publication[order.dendrogram(euc_dend)]
labels(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.") 
dend_cols <- paste0(sample_info_tab$case, sample_info_tab$publication)[order.dendrogram(euc_dend)]
labels(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")
#yikes, Michaolich is still batch effecting hard
#ASK MIKE LEE ABOUT FIXING THIS

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab2, taxa_are_rows=T)
vst_tax_phy <- tax_table(bac_pseq)
vst_physeq <- phyloseq(vst_count_phy, vst_tax_phy,sample_data(bac_pseq))


vst_physeq
# # generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
# 
library(ggplot2)
plot_ordination(vst_physeq, vst_pcoa, color="publication") + 
   geom_point(size=1) + labs(col="case") + 
   coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") 

#####ok lets run some stats on the normalized counts##########
library(tidyverse)
library(matrixTests)
case_sum<-psmelt(vst_physeq)%>%
  select(Sample,OTU,case,Abundance)%>%
  group_by(OTU,case)%>%
  pivot_wider(id_cols = c(Sample,case), 
              names_from = OTU,
              values_from = Abundance)

dim(case_sum)
mat<-case_sum[,3:length(case_sum)]
krus_case<-col_kruskalwallis(x=mat,g = case_sum$case)
welch_case<-col_oneway_welch(x=mat, g=case_sum$case)
krus_case_sig<-krus_case%>%filter(pvalue<0.05)%>%arrange(pvalue)
welch_case_sig<-welch_case%>%filter(pvalue<0.05)%>%arrange(pvalue)
colnames(krus_case_sig)
colnames(welch_case_sig)
dim(krus_case_sig)
dim(welch_case_sig)
write.table(krus_case,"krus_case.tsv",sep = "\t")
write.table(welch_case,"welch_case.tsv",sep = "\t")
sig_names<-intersect(rownames(krus_case_sig),rownames(welch_case_sig))
length(sig_names)
sig_names
#########################################################
#ok now lets use our findings from the stat test in order to
#filter out non significant case assocaited taxa
bac_pseq_prune <- prune_taxa(taxa = sig_names, bac_pseq_no_neg_unique)
# bac_pseq_prune # [ 1353 taxa and 141 samples ]
# bac_pseq_prune <- prune_samples(sample_sums(bac_pseq_prune) > 1, bac_pseq_prune)
# bac_pseq_prune <- prune_taxa(taxa_sums(bac_pseq_prune) > 1, bac_pseq_prune)
bac_pseq_prune # [ 2089 taxa and 126 samples ]:


#filter out the suppper high stuff like molecular function and biological process
# [1] "GO:0045212-obsolete neurotransmitter receptor biosynthetic process"                                                     
# [2] "GO:0001319-obsolete inheritance of oxidatively modified proteins involved in replicative cell aging"                    
# [3] "GO:1900008-obsolete negative regulation of extrachromosomal rDNA circle accumulation involved in replicative cell aging"
# [4] "GO:0008150-biological_process"                                                                                          
# [5] "GO:0001302-obsolete replicative cell aging"                                                                             
# [6] "GO:0003674-molecular_function"
#NAILLLLED IT

#GOTERM DEPTH LEVEL1
# [1] "GO:0051704-multi-organism process"                     "GO:0008152-metabolic process"                         
# [3] "GO:0051703-intraspecies interaction between organisms" "GO:0110148-biomineralization"                         
# [5] "GO:0043473-pigmentation"                               "GO:0002376-immune system process"                     
# [7] "GO:0044419-interspecies interaction between organisms" "GO:0048511-rhythmic process"                          
# [9] "GO:0023052-signaling"                                  "GO:0007610-behavior"                                  
# [11] "GO:0015976-carbon utilization"                         "GO:0051179-localization"                              
# [13] "GO:0040011-locomotion"                                 "GO:0032502-developmental process"                     
# [15] "GO:0032501-multicellular organismal process"           "GO:0040007-growth"                                    
# [17] "GO:0009758-carbohydrate utilization"                   "GO:0022610-biological adhesion"                       
# [19] "GO:0022414-reproductive process"                       "GO:0050896-response to stimulus"                      
# [21] "GO:0006791-sulfur utilization"                         "GO:0098754-detoxification"                            
# [23] "GO:0000003-reproduction"                               "GO:0065007-biological regulation"                     
# [25] "GO:0009987-cellular process"                           "GO:0019740-nitrogen utilization"                      
# [27] "GO:0032947-molecular adaptor activity"                 "GO:0060090-translation regulator activity"            
# [29] "GO:0045182-protein folding chaperone"                  "GO:0044183-structural molecule activity"              
# [31] "GO:0005198-molecular carrier activity"                 "GO:0140104-catalytic activity"                        
# [33] "GO:0003824-toxin activity"                             "GO:0090729-molecular transducer activity"             
# [35] "GO:0060089-cargo receptor activity"                    "GO:0038024-molecular function regulator"              
# [37] "GO:0098772-transporter activity"                       "GO:0005215-small molecule sensor activity"            
# [39] "GO:0140299-protein tag"                                "GO:0031386-antioxidant activity"                      
# [41] "GO:0016209-binding"                                    "GO:0005488-nutrient reservoir activity"               
# [43] "GO:0045735-multi-organism process


bac_pseq_prune<- subset_taxa(as.numeric(data.frame(tax_table(bac_pseq_prune))$depth)>1, physeq = bac_pseq_prune) 
bac_pseq_prune <- prune_samples(sample_sums(bac_pseq_prune) > 10, bac_pseq_prune)
bac_pseq_prune <- prune_taxa(taxa_sums(bac_pseq_prune) > 10, bac_pseq_prune)
bac_pseq_prune#[ 2077 taxa and 126 samples ]:
# bac_pseq_prune_deeper<- prune_taxa(taxa_sums(bac_pseq_prune_deeper) > 10, bac_pseq_prune_deeper)
# bac_pseq_prune_deeper<- prune_samples(sample_sums(bac_pseq_prune_deeper) > 10, bac_pseq_prune_deeper)
# bac_pseq_prune_deeper # otu_table()   OTU Table:         [ 50 taxa and 71 samples ] not bad
#THIS IS A HACK AND i NEED TO FIND A BETTER WAY TO VECTORIZE THE NAMES INSTEAD OF AGGLOMERATING THE GO TERMS BY NAMES AND LOSING THE RESOLUTION
#OF HAVING BEING ABLE TO SWAP THE GO TERMS WITH THE NAMES IN THE VISUALIZATION WITHOUT LOSING THE RESOLUTION OF AGGLOMERATING GO TERMS WITH THE SAME NAME
#AGLOMERATE BY NAME BECAUSE THERRE ARE ONLY 1 GO TAG THAT HAVE THE SAME NAME
#bac_pseq_prune<-tax_glom(bac_pseq_prune, taxrank = "name")
#CHANGE THE GO TERM NAME FOR THE IDENTIFIER

#UPDATE
#I fixed the tax_aglom naming issue but appending the goterm identifier to the name, making everything unique again
names<-paste(taxa_names(bac_pseq_prune),get_taxa_unique(bac_pseq_prune,taxonomic.rank = "name" ),sep = "-")
taxa_names(bac_pseq_prune)<-names#paste0(taxa_names(bac_pseq_prune),get_taxa_unique(bac_pseq_prune,taxonomic.rank = "name" )

# bac_pseq_prune2<- subset_taxa(as.numeric(data.frame(tax_table(bac_pseq_no_neg_unique))$depth)==1, physeq = bac_pseq_no_neg_unique) 
# bac_pseq_prune2
# names<-paste(taxa_names(bac_pseq_prune2),get_taxa_unique(bac_pseq_prune2,taxonomic.rank = "name" ),sep = "-")
# taxa_names(bac_pseq_prune2)<-names
# taxa_names(bac_pseq_prune2)

####################
#DMM MoDELING TIME
####################

library(microbiome)
library(DirichletMultinomial)
library(phyloseq)
#convert counts to a matrix
dat <- abundances(bac_pseq_prune)
count <- as.matrix(t(dat))

lvls <- c("Control_Healthy", "Control_Sick", "COVID19")
pheno<-factor(sample_data(bac_pseq_prune)$case, levels=lvls)
pheno


length(taxa_names(bac_pseq_prune))

#Fit the model

bestgrp <- dmngroup(count, pheno, k=1:5, verbose=TRUE,simplify = TRUE, mc.preschedule=FALSE)

#fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)
fit<-dmngroup(dmn,count, k=1:5, simplify = TRUE,
           .lapply = parallel::mclapply)
#fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)

#Check the model fit with different number of mixture componenets using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

## best fit for groups 'Lean' and 'Obese'; full example in vignette.
bestgrp <- dmngroup(count, pheno, k=1:8, verbose=TRUE, 
mc.preschedule=FALSE)


data(fit)
## End(Not run)
data(bestgrp)
bestgrp
bestgrp[["Obese"]]

#identify the number of clusters that best fits the model
best <- fit[[which.min(lplc)]]
best
best <-fit[[3]]
#save.image(file = "bac_go_terms_dmm.rdata")

#make a heatmap visualization of the cluster
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =log2, lblwidth = 0.2 * nrow(count))
#make a heatmap visualization of the cluster
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =sqrt, lblwidth = 0.2 * nrow(count))
#print out the theta values
mixturewt(best)
#save datasheet that show which GO terms contributed to each dmm group
write.table(fitted(best),"combined_bac_GO_TERMS_DMM_contributions.tsv", sep="\t")
#save a datasheet that identifies which sample belongs to which dmm group
ass <- apply(mixture(best), 1, which.max)
write.table(ass,"combined_bac_GO_TERMS_DMM_groups.tsv",sep="")

#make a copy of the phyloseq object so you dont jack it up
physeq<-bac_pseq_prune
meta<-data.frame(sample_data(physeq))
meta$dmm<-ass
library(ggplot2)
library(dplyr)
library(reshape2)
#Go through each dmm cluster and display Go_term contributions to each cluster
for (k in seq(ncol(fitted(best))))
{
  d <- melt(fitted(best))
  colnames(d) <- c("GO", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange GOs by assignment strength
    arrange(value) %>%
    mutate(GO = factor(GO, levels = unique(GO))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.99))
  
  p <- ggplot(d, aes(x = GO, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers in  : GO  Terms cluster type ", k, sep = ""))+
           ggsave(filename = paste("contribution_bar_plot_",k,sep = ""),device = "png")
  #paste(p,k, sep = "")<-p
  #print(k)
  print(p)
  #print(paste(p,k, sep=""))
}


#add the dmm group to the metadata
#sample_data(bac_pseq_prune_deeper2)$dmn<-ass
# #melt the phyloseq object into tidy form
#  tmp<-psmelt(bac_pseq_prune_deeper2)
#  tmp<-as_tibble(tmp)
# 
# #subset the dataset to only include the case, Go_term, count, and dmm group.
# #obtain the avergage count for each Go term
# #order the go terms from hight to lowest count
# 
# d2<-tmp %>%
#   select(case,name,Abundance, dmn)%>%
#   group_by(name,case, dmn) %>%
#   summarise(avg = mean(Abundance)) %>%
#   arrange(desc(avg))
# 
save.image(file = "combined_bac_go_terms_dmm.rdata")
# 
# #move each dmm group into a colum of its own
# d3<-tidyr::spread(d2,dmn, avg)
# 
# #get the total count of the go terms and oder from greates to lowest
# d3[3:5]
# d3$tot<-rowSums(d3[3:5], na.rm = T)
# d3<-d3%>%arrange(desc(tot))
# d3$tot<-NULL
# 
# d3<-d3%>%gather(data = d3,avg,3:5)
# colnames(d3)<-c("name","case", "dmn","avg")
# d4<-d3%>%top_n(100, c(avg))
# unique(d4$name)
# library(ggpubr)
# ggballoonplot(d4, y ="name",x = "case", size = "avg", facet.by = "dmn",fill = "avg")+scale_fill_viridis_c()
# 
#bac_pseq_prune_deep
#add the dmm group to the metadata
sample_data(bac_pseq_prune)$dmn<-ass
tmp<-psmelt(bac_pseq_prune)
tmp<-as_tibble(tmp)

d2<-tmp %>%
  select(case,name,Abundance, dmn)%>%
  group_by(name,case, dmn) %>%
  summarise(avg = mean(Abundance),std = sd(Abundance)) %>%
  group_by(name, case,dmn, std)%>%
  arrange(desc(std))
d2
#d3<-tidyr::spread(d2,dmn, avg)
d3<-tidyr::spread(d2,case,avg)
library(tidyverse)
d3[4:6]

#library(mosaic)
#d3$tot<-sd(d3[3:5])
#d3$tot<-rowSums(d3[3:5], na.rm = T)
d3<-d3%>%arrange(desc(std))
rank<-as.character(d3$name)
dim(d3)
d4<-d3[1:100,]

d4<-d4%>%
  gather(data = d4,avg,4:6)

colnames(d4)<-c("name","dmn","std", "group","avg")

#d4<-d3%>%group_by(name,group)
d4<-d4%>%filter(!is.na(avg))
library(ggpubr)
d4$dmn<-as.character(d4$dmn)
ggballoonplot(data = d4, y ="name",x = "dmn", size = "avg", facet.by = "group",fill = "avg")+scale_fill_viridis_c()
write.table(meta(bac_pseq_prune),file = "combined_pseq_meta.tsv",sep = '\t')
library(mosaic)


my_tbl<-tally(~dmn,data =meta(bac_pseq_prune),format = 'count')
my_tbl
library(microbiome)
my_tbl<-tally(~case+ publication+dmn,data =meta(bac_pseq_prune),format = 'count')
my_tbl
a<-as_tibble(my_tbl)

ggballoonplot(data = a, y ="case",x = "dmn", facet.by = "publication",size = "n", fill = "n")+
  scale_fill_viridis_c(option = "C")+
  ggsave(filename = "pub_vs_disease_vs_dmm_dotplot.png",
         device = "png",
         #width = "8",
         #height = "6",
         #units = "in",
         dpi = 600,
         path = "D:/github/microbial/Rdata/GO_term_analysis/Figures/")

my_tbl<-tally(~case+ dmn,data =meta(bac_pseq_prune),format = 'count')
my_tbl_prop<-tally(~case+ dmn,data =meta(bac_pseq_prune),format = 'prop')
my_tbl
my_tbl_prop<-as_tibble(my_tbl_prop)
b<-aov(n~dmn*case,my_tbl_prop)
msummary(a)
c<-TukeyHSD(b)
chisq.test(my_tbl)
#mplot(b$model)
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="case",x = "dmn", size = "n", fill = "n")+
  scale_fill_viridis_c(option = "B")+
  ggsave(filename = "disease_vs_dmm_dotplot.png",
         device = "png",
         #width = "8",
         #height = "6",
         #units = "in",
         dpi = 600,
         path = "D:/github/microbial/Rdata/GO_term_analysis/Figures/")

my_tbl<-tally(publication ~ dmn,data =meta(bac_pseq_prune),format = 'count')
#my_tbl<-tally(publication ~ dmn,data =meta(bac_pseq_prune),format = 'prop')
a<-as_tibble(my_tbl)
b<-aov(n~publication,a)
msummary(b)
msummary(c)
c<-TukeyHSD(b)
c
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="publication",x = "dmn", size = "n", fill = "n")+scale_fill_viridis_b(option = "B")
my_tbl<-tally(publication ~ dmn+case,data =meta(bac_pseq_prune),format = 'count')
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="publication",facet.by = "case",x = "dmn", size = "n", fill = "n")+scale_fill_viridis_b(option = "B")

###########################################
###Dont forget to save you shit HERE#######
###########################################
#save.image(file = "bac_go_terms_dmm.rdata")
#

library(dplyr)
count<-abundances(bac_pseq_prune)
select <- order(rowMeans(count),decreasing=TRUE)
select2<-log1p((count)[select,])
tmp<-rownames(select2)
dim(select2)

#select2$mol<-tmp
select3<-as.table(select2)

#library(matrixTests)
#library(genefilter)
#select3<-as_tibble(select2)%>%summarise(std=rowFtests(select2))%>% arrange(desc(std))
#select2$mol<-tmp
select2<-as_tibble(select2)
select2
#select2$mean<-rowMeans(select2)
rownames(select2)<-tmp
rownames(select2)
sam<-data.frame(sample_data(bac_pseq_prune))
df<-as.data.frame(sample_data(bac_pseq_prune))
df<-as_tibble(df)
df<-df%>%select(dmn, publication, sample_type,case)#,dmn,body_site)
df
df<-as.data.frame(df)
df

#select2<-sqrt((count)[select,])
dim(select2)
colnames(select2) <- colnames(otu_table(bac_pseq_prune))
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
df
colnames(df)<-c("dmm_cluster", "Publication","Sample_Type","Case")
# Specify colors
ann_colors = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF"),#,"4"="#631879FF","5"="#008280FF"),
  Publication=c("Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Sample_Type =c("COVID_19"="firebrick",
                 "Healthy"="forestgreen",
                 "Community_acquired_pneumonia"="dodgerblue4",
                 "Obese_Asthma"="goldenrod3",
                 "Obese_Smoker"="goldenrod4", 
                 "Obese"="goldenrod1",
                 "Asthma"= "darkorange2", 
                 "Asthma_Smoker"="darkorange4",
                 "Asthma_Ex_smoker"="darkorange3",
                 "Smoker"="gray27",
                 "Obese_Asthma_Smoker"="black"),
  Case=c("COVID19"="firebrick", "Control_Healthy"="forestgreen","Control_Sick"="dodgerblue4"))
rowMeans(select2)

select3<-select2%>%filter(rowMeans(select2)>0.1)
xx <- pheatmap(mat = select3,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),
               annotation_col=df,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
              annotation_row = df_row)
xx <- pheatmap(mat = select2,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10),
               annotation_col=df,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               annotation_row = df_row,
               cluster_row = F,
               cluster_col=T)
