#Data Wrangling
#COVIRT19 microbial subgroup seqscreen analysis
#Date : 29 MARCH 2021
#Maintained by :
#jochum, Michael D. 
#Baylor College of Medicine 
#michael.jochum@bcm.edu


###########################################################################################
###########################PRE-PROCESSING##################################################
###########################################################################################
# The purpose of this code is to:
#   take the raw seqscreen GO Term counts and convert them into a working phyloseq object
# conduct some preprocessing on the phyloseq object that:
#   filters out batch effect samples and GO terms with little to no abundance

#load the libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
setwd("C:/github/microbial/GO_Terms_scripts/")

#import the count and metadata tables 
raw<-as_tibble(read.table("datasets/Combined_BALF_GO_Terms_parent_propagated_03292021.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
meta<-as.data.frame(read.table("datasets/Combined_BALF_GO_Terms_metadata_03292021.tsv",header = T, sep = "\t",row.names = 1))


#select the count columns
#pivot them by sample,type, and abundance
#filter Terms withzero counts values 
#pivot the samples back into the rows, putting a 0 in place of NA values

df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)

#convert the depth to a character class
df$depth<-as.character(df$depth)

term<-df%>%filter(type!="NA")%>%filter(type%in%c("bac","arc","vir"))%>%group_by(GO_term,namespace,depth,name)%>%
  summarise(across(.cols = where(is.numeric), sum))

term_tax<-term%>%select(GO_term,namespace,depth,name)
term_tax<-data.frame(term_tax, row.names=1)
term_counts<-data.frame(term[5:172], row.names = term$GO_term)

term_counts_phy <- otu_table(term_counts, taxa_are_rows=TRUE)
term_tax_phy <- tax_table(as.matrix(term_tax), errorIfNULL=TRUE)
term_sam<-sample_data(meta)

#make the phyloseq object
term_pseq <- phyloseq(term_counts_phy, term_tax_phy, term_sam)
term_pseq# [ 14581 taxa and 167 samples ] [ 27077 taxa and 167 samples ]

#remove terms Biological process and Molecular Function depth zero
filtme<-c("GO:0003674","GO:0008150") 
term_pseq <- prune_taxa(taxa=!taxa_names(term_pseq)%in%filtme, term_pseq)
term_pseq #[ 14579 taxa and 167 samples ]



####################################################
####Sample removal stage 
####################################################
#remove the negative control samples
term_pseq_no_neg<-subset_samples(term_pseq, sample_type!="neg_control")
term_pseq_no_neg# [ 14579 taxa and 162 samples ]:
#remove the samples with Unknown and Control Sick cases because they dont provide us with anything of value
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, sample_type!="Unknown")
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, case!="Control_Sick")
term_pseq_no_neg# [ 14597 taxa and 105 samples ]
# remove the remaining Michalovich that are producing an strong batch effect
term_pseq_no_neg<-subset_samples(term_pseq_no_neg,publication!="Michalovich")
term_pseq_no_neg# [ 14597 taxa and 102 samples ]
# remove the samples that were viral enriched and also producing an strong batch effect
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, bioproject!="PRJNA605907")
#term_pseq_no_neg# [ 14597 taxa and 86 samples ]
#remove GO Terms with zero counts
term_pseq_no_neg<-prune_taxa(taxa = taxa_sums(term_pseq_no_neg)>0,x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ]
#remove samples with zero counts
term_pseq_no_neg<-prune_samples(samples = sample_sums(term_pseq_no_neg)>0,x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ] # [ 25426 taxa and 86 samples ]
#make a copy of the phyloseq_object with the GO_TERM TAGS
term_pseq_no_neg_gonames<-term_pseq_no_neg

#GO_term filtering
#Remove GO TERMS with less than 1% prevalance and 1% abundance
term_pseq_no_neg_core<-core(x = term_pseq_no_neg,detection = 1/100,prevalence = 1/100)



####################################################
####TAKE THIS OUT BC MIKE SAID TO MAKE THE TAGS CLEAN
####################################################
# #make a copy of the phyloseq_object with the GO_TAGS
# term_pseq_no_neg_gonames<-term_pseq_no_neg
# 
# #rename the taxa by the GO_term descriptions
# tax<-data.frame(tax_table(term_pseq_no_neg))
# names<-paste(rownames(tax),tax$name,sep="-")
# taxa_names(term_pseq_no_neg)<-names


###########################################################################################
###########################BETADISPER##################################################
##############Multivariate homogeneity of groups dispersions (variances)#############
###############################################################################################
#######################################BETADISPER#####################################

#Cases

#lets check for homogeneity and within group sample distances to make sure we dont have anymore stong batch effects
library(DESeq2)
deseq_counts <- phyloseq_to_deseq2(term_pseq_no_neg_core,design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
dis <- dist(t(vst_trans_count_tab))
sam<-meta(term_pseq_no_neg_core)

library(vegan)

## Calculate multivariate dispersions using beatdisper testing
mod <- betadisper(dis, sam$case)
mod
## Perform an Analysis of Variance Table test on the betadisper
summary(anova(mod))
## Perform a PERMANOVA Permutation test for F
permutest(mod, pairwise = FALSE, permutations = 999)
## Perform Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
mod.HSD

vec<-as.data.frame(mod$vectors[,1:2])
vec$case<-sam$case
library(ggsci)
library(ggpubr)

ggscatter(data = vec,
          x = "PCoA1",
          y = "PCoA2",
          color = "case",
          title = "PCoA of within group homogeneity (Case)",
          palette = "npg",
          size = 1,
          rug = T,
          conf.int = T,
          ggtheme = theme_bw(),ellipse.border.remove = F,
          ellipse = T,ellipse.type = "confidence" ,ellipse.alpha = 0,
          cor.method = "spearman",
          fullrange = T,
          star.plot = T)+
  coord_fixed(sqrt(mod$eig[2]/mod$eig[1]))
  
## Draw a boxplot of the distances to centroid for each group
boxplot(mod)
# # making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sam)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
# 
# # generating and visualizing the PCoA with phyloseq
 vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
 eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
eigen_vals
 # 
plot_ordination(vst_physeq, vst_pcoa, color="patient") + 
 geom_point(size=1) + labs(col="case") + 
#   geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
 coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_npg()+theme_bw()



deseq_counts <- phyloseq_to_deseq2(term_pseq_no_neg_core,design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
dis <- dist(t(vst_trans_count_tab))
sam<-meta(term_pseq_no_neg_core)
## Calculate multivariate dispersions using beatdisper testing
mod <- betadisper(dis, sam$patient)
mod
## Perform an Analysis of Variance Table test on the betadisper
anova(mod)
## Perform a PERMANOVA Permutation test for F
permutest(mod, pairwise = FALSE, permutations = 999)
## Perform Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

vec$patient<-sam$patient
vec$case<-sam$case
library(ggsci)
library(ggpubr)

ggscatter(data = vec,
          x = "PCoA1",
          y = "PCoA2",
          color = "patient",shape="case",
          title = "PCoA of within group homogeneity (patient)",
#          palette = "npg",
          size = 1,
          rug = T,
          conf.int = T,
          ggtheme = theme_bw(),ellipse.border.remove = F,
          ellipse = T,ellipse.type = "confidence" ,ellipse.alpha = 0,
          cor.method = "spearman",
          fullrange = T,
          star.plot = T)+
  coord_fixed(sqrt(mod$eig[2]/mod$eig[1]))

## Draw a boxplot of the distances to centroid for each group
boxplot(mod)
# # making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sam)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
# 
# # generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
eigen_vals
# 
plot_ordination(vst_physeq, vst_pcoa, color="patient") + 
  geom_point(size=1) + labs(col="case") + 
  #   geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_npg()+theme_bw()
 


# 

#Save the phyloseq object
#save.image(file = "./images/0_preprocessing.RDA")

###########################PRE-PROCESSING ENDS HERE#######################################



###########################################################################################
########################### MAAslin2 ##################################################
###########################################################################################
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(microbiome)
library(parallel)
setwd("C:/github/microbial/GO_Terms_scripts/")

#intermediate analysis of the shen samples
#library(ggpubr)
# a<-as.data.frame(sample_sums(term_pseq_no_neg))
# a$id<-rownames(a)
# a<-as_tibble(a)%>%arrange(desc(value))
# b<-meta(term_pseq_no_neg)
# b$id<-rownames(b)
# 
# a
# b
# c<-full_join(a,b)
# c$sample_sums<-c$`sample_sums(term_pseq_no_neg)`
# c$`sample_sums(term_pseq_no_neg)`<-NULL
# shen<-c%>%filter(publication=="Shen",case=="COVID19")%>%select(patient,sample_sums,id)%>%arrange(desc(sample_sums))
# gghistogram(data = shen,
#             x = "sample_sums",
#             rug = T,
#             add_density = F,
#             label = "id",
#             fill = "patient",
#             alpha=0.2,
#             palette = "aaas")+
#   xscale(.scale = "log10",.format = T)+
#   facet_wrap(nrow =2,ncol = 4,facets = ~patient)
# write.table(shen,"shen_sample_sums.tsv",sep = "\t",quote = F,row.names = F)

<<<<<<< Updated upstream

=======
Sys.Date()
>>>>>>> Stashed changes


#Transform GO Terms counts to relative abundances
term_pseq_no_neg_comp<-microbiome::transform(x = term_pseq_no_neg,transform = "compositional")
term_pseq_no_neg_comp#[ 13534 taxa and 86 samples ]
#lets try filtering terms upstream using a 1% abundance and 1% prevalence cutoffs
term_pseq_no_neg_comp_core<-core(x = term_pseq_no_neg_comp,detection = 1/100,prevalence = 1/100)
term_pseq_no_neg_comp_core #[ 90 taxa and 86 samples ]
term_pseq_no_neg_comp_core #[ 64 taxa and 102 samples ] with viral enriched shen samples

#convert the metadata compositionally transformed count table into a dataframe
df_input_data<-data.frame(t(otu_table(term_pseq_no_neg_comp_core)))
df_input_metadata<-data.frame(sample_data(term_pseq_no_neg_comp_core))

#Multivariable Association with Linear Models
#run a three way comparison of cases using COVID19 as the reference group 
#Control for the random effects observed from using multiple samples from the sample publication and patient

#use a minimum of 1% abundance and 10% prevalence cutoff 
#EDIT:heads up I removed the minimum abundance and prevlance cutoffs in this analysis because I did them upstream


#filter Terms that do not have a p value less that 0.05 using benjamini-hochberg multiple test correction
case<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./results/terms_vs_case_comp",
  min_abundance = 0, 
  min_prevalence = 0, 
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("patient","publication"),
  fixed_effects = c("case"),
  correction="BH",
  standardize = TRUE,
  cores = detectCores(),
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)


case_res<-case$results%>%filter(pval<=0.05)
length(unique(case_res$feature)) #66 before #61 with including the viral enriched samples
case_res

library(mosaic)
term_pseq_no_neg_survival<-subset_samples(physeq = term_pseq_no_neg,survival!="NA")
term_pseq_no_neg_survival #[ 13534 taxa and 25 samples ]
#Transform GO Terms counts to relative abundances
term_pseq_survival_no_neg_comp<-microbiome::transform(x = term_pseq_no_neg_survival,transform = "compositional")
term_pseq_survival_no_neg_comp#[ 13534 taxa and 25 samples ]
#lets try filtering terms upstream using a 1% abundance and 1% prevalence cutoffs
term_pseq_survival_no_neg_comp_core<-core(x = term_pseq_survival_no_neg_comp,detection = 1/1000,prevalence = 1/1000)
term_pseq_survival_no_neg_comp_core # [ 435 taxa and 25 samples ]


#convert the metadata compositionally transformed count table into a dataframe
COVID19_df_input_data<-data.frame(otu_table(term_pseq_survival_no_neg_comp_core))
COVID19_df_input_metadata<-data.frame(sample_data(term_pseq_survival_no_neg_comp_core))

# 
# library(dplyr)
# library(mosaic)
# tally(~accession+survival+patient+publication,COVID19_df_input_metadata,"data.frame")%>%filter(Freq>0)
# #heads up I removed the minimum abundance and prevlance cutoffs in this analysis because I did them upstream
survival<-Maaslin2(
  input_data = COVID19_df_input_data,
  input_metadata = COVID19_df_input_metadata,
  output="./results/survival",
  min_abundance = 0,
  min_prevalence = 0,  #test only features with at least 10% non-zero values. min 2.5samples to be run
  normalization = "CLR",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.25,
  random_effects = c("patient"),
  fixed_effects = c("survival"),
  correction="BH",
  standardize = TRUE,
  cores = detectCores(),
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =num)


survival_res<-survival$results%>%filter(pval<=0.05)
length(unique(survival_res$feature)) 
taxa_names(term_pseq_no_neg_comp_core)

