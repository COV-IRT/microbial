library(tidyverse)
library(phyloseq)
getwd()
#setwd("c:/github/microbial/step2_kraken2_analysis/")

# if wanting to load the objects
#load("heat-tree.RData")

# reading in metadata
sample_info_tab <- read.table("Combined_BALF_GO_Terms_metadata.tsv", header = TRUE, sep = "\t")
dim(sample_info_tab)
head(sample_info_tab)[, 1:5]

# reading in tax file
full_tax_tab <- read.table("all_step2_kraken2_summaries.tsv", header = TRUE, row.names = 1, sep = "\t", quote = '', comment.char="")
dim(full_tax_tab)
names(full_tax_tab)

# making table of just counts (taxid as rows)
tax_counts_tab <- full_tax_tab[, grep("counts", names(full_tax_tab))]
head(tax_counts_tab)[, 1:5]
# shortening names
names(tax_counts_tab) <- sub("_step2_kraken2_mason_db_read_counts", "", names(tax_counts_tab))
head(tax_counts_tab)[, 1:5]

# making table of just tax info
tax_info_tab <- full_tax_tab[, 1:7]
head(tax_info_tab)

# filtering down to only the 86 samples in the GO term analysis
target_samples <- scan("GO_TERM_samples.tsv", what = "character")[-1]
target_samples_no_neg<-target_samples
#add the negative controls back into the analysis
negs<-c("CRR125995",'CRR125996',"CRR125997","CRR125998")
target_samples<-c(target_samples,negs)
head(target_samples)
length(target_samples)

target_tax_counts_tab <- tax_counts_tab %>% select(all_of(target_samples))
dim(target_tax_counts_tab)

# another look they were all found even though right size table came through
setdiff(target_samples, names(target_tax_counts_tab))

target_sample_info_tab <- sample_info_tab %>% filter(accession %in% target_samples)
dim(target_sample_info_tab)
setdiff(target_samples, target_sample_info_tab$accession)

target_sample_info_tab <- target_sample_info_tab %>% column_to_rownames("accession")
# filtering down to just Bacteria
bac_tax_counts_tab <- target_tax_counts_tab[!is.na(tax_info_tab$domain) & tax_info_tab$domain == "Bacteria", ]
dim(bac_tax_counts_tab)
colSums(bac_tax_counts_tab) / colSums(target_tax_counts_tab) * 100
# that is removing a ton, like 99% from some, i hope it's just that there are a ton of unclassifieds like kmer tax-assignment does, let's check though

target_tax_counts_tab[1, ] / colSums(target_tax_counts_tab) * 100
# oh geez, yea
round(target_tax_counts_tab[1, ] / colSums(target_tax_counts_tab) * 100, 3) %>% as.matrix %>% as.vector %>% summary
# over 75% of them are >88% Unclassified at all, silly short-read taxonomy, ha

colSums(bac_tax_counts_tab) / (colSums(target_tax_counts_tab) - target_tax_counts_tab[1, ]) * 100
# there are still some with very low relative proportions, like CRR125941 CRR125949 only having 7% remaining, let's see what's in them
target_tax_counts_tab %>% rownames_to_column("taxid") %>% select(taxid, CRR125941, CRR125949) %>% arrange(desc(CRR125949)) %>% head

tax_info_tab[row.names(tax_info_tab) == "2697049", ] # ahh, a virus, are these viral enriched samples?

# moving onto the figure for now
bac_tax_info_tab <- tax_info_tab[!is.na(tax_info_tab$domain) & tax_info_tab$domain == "Bacteria", ]
dim(bac_tax_info_tab)
# making sure all is good
all.equal(row.names(bac_tax_info_tab), row.names(bac_tax_counts_tab))



#MJ-Make a pseq-object so you can have something to filter in and out of for decontam and metacoder...sorry Mike =)
pseq<-phyloseq(otu_table(bac_tax_counts_tab,taxa_are_rows = T),tax_table(as.matrix(bac_tax_info_tab)),sample_data(target_sample_info_tab))
pseq

library(decontam)

#export the shen et al. samples for parsing with decontam
pseq_shen<-subset_samples(pseq,bioproject=="CRA002476")
pseq_not_shen<-subset_samples(pseq,bioproject!="CRA002476")
pseq_shen

count_df_shen<-abundances(pseq_shen)
tax_table(pseq_shen) %>% as("matrix")# %>% subset(select = colSums(!is.na(.)) > 0) %>% as_tibble(rownames = "OTU") -> taxon.tbl
library(microbiome)
tax_df_shen<-data.frame(tax_table(pseq_shen))
tax_df_shen
sam_shen<-meta(pseq_shen)
vector_for_decontam <- sam_shen$case=="Control_Neg"
vector_for_decontam                    

colSums(count_df_shen)
count_df_shen<-as.data.frame(otu_table(pseq_shen))
colSums(count_df_shen)

contam_df_shen <- isContaminant(t(count_df_shen), neg=vector_for_decontam)
table(contam_df_shen$contaminant) # identified 4644 as contaminants 
# getting vector holding the identified contaminant IDs
contam_shen <- row.names(contam_df_shen[contam_df_shen$contaminant == TRUE, ])
contam_shen<-as.factor(contam_shen)
contam_shen



count_df_shen<-as_tibble(count_df_shen,rownames = "taxid")
count_df_shen
#count_df_shen[1]<-NULL

clean_count_df_shen<-count_df_shen%>%filter(!taxid%in%contam_shen)

a<-clean_count_df_shen%>%select(-taxid)
b<-count_df_shen%>%select(-taxid)
perc_remain<-(colSums(a)/colSums(b))*100
#write.table(perc_remain,"percent_remaining_post_decontam.tsv",sep = "\t",quote = F,row.names = F)
tax<-as_tibble(bac_tax_info_tab,rownames="taxid")
contaminants<-tax%>%filter(taxid %in% contam_shen)
not_contaminants<-tax%>%filter(!taxid %in% contam_shen)

contaminants%>%filter(genus=="Variovorax")
not_contaminants%>%filter(genus=="Variovorax")



not_contaminants
#write.table(x = contaminants,file = "Shen_contaminants.tsv",sep = "\t",row.names = F,quote = F)
#write.table(x = not_contaminants,file = "Shen_not_contaminants.tsv",sep = "\t",row.names = F,quote = F)

contam_shen<-as.factor(contam_shen)

contaminants$taxid<-as.character(contaminants$taxid)
not_contaminants$taxid<-as.character(not_contaminants$taxid)

pseq_shen_decontam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)
pseq_shen_decontam


pseq_decontam<-merge_phyloseq(pseq_not_shen,pseq_shen_decontam)
pseq_decontam_no_neg<-subset_samples(physeq = pseq_decontam,case!="Control_Neg")

pseq_decontam_no_neg_core<-core(x = pseq_decontam_no_neg,detection = 0,prevalence = 50/100)
pseq_decontam_no_neg_core #[ 566 taxa and 86 samples ]






#install.packages('metacoder')
library(metacoder)
# making our data look like theirs in this example page: https://grunwaldlab.github.io/metacoder_documentation/example.html
  # which we can look at in the hmp_otus object after loading the metacoder library
hmp_otus
# they have an ID column, a lineage column, then the counts, we're going to turn our tax info into their lineage column format
hmp_otus$lineage %>% head(30) %>% tail

lineage_vector <- paste0(
    "d__", bac_tax_info_tab$domain, ";",
    "p__", bac_tax_info_tab$phylum, ";",
    "c__", bac_tax_info_tab$class, ";",
    "o__", bac_tax_info_tab$order, ";",
    "f__", bac_tax_info_tab$family, ";",
    "g__", bac_tax_info_tab$genus, ";",
    "s__", bac_tax_info_tab$species)

# and the don't have the "NA" ones, they just have nothing for them, so removing that business so ours stop at the last rank present like theirs
lineage_vector <- gsub("(;[c-s]__NA)*$", "", lineage_vector)
head(lineage_vector)
# beautiful, building table like hmp_otus
tab_for_taxa_object <- data.frame(taxid = row.names(bac_tax_info_tab), lineage = lineage_vector, bac_tax_counts_tab)
dim(tab_for_taxa_object)
tab_for_taxa_object<-tab_for_taxa_object%>%select(-negs)

head(tab_for_taxa_object)
head(hmp_otus)[, 1:10]
str(tab_for_taxa_object)
# i think all's ok, let's give it a shot
sample_data(pseq_decontam_no_neg_core)$case<-gsub('Control_Healthy','Uninfected',sample_data(pseq_decontam_no_neg_core)$case)
sample_data(pseq_decontam_no_neg_core)$case<-gsub('COVID 19','COVID-19',sample_data(pseq_decontam_no_neg_core)$case)
sample_data(pseq_decontam_no_neg_core)$case<-gsub('Control_Sick','CAP',sample_data(pseq_decontam_no_neg_core)$case)







pseq_obj<-parse_phyloseq(pseq_decontam_no_neg_core)

taxa_obj <- parse_tax_data(tab_for_taxa_object,
                           class_cols = "lineage",
                           class_sep = ";",
                           class_regex = "^([c-s]{1})__(.*)$",
                           class_key = c(tax_rank = "info", tax_name = "taxon_name"))


library(taxa)
# looks ok so far, following along with same example page to filter low abundance things, just picking anything
taxa_obj_filt <- filter_taxa(taxa_obj, n_obs > 5)
pseq_obj_filt <- filter_taxa(pseq_obj, n_obs > 5)


# this is where we normalize across samples (i'm curious what this spits out when there are negative values in the table like from the vst, ha)
#################################################################################
pseq_obj_filt$data$tax_proportions <- calc_obs_props(pseq_obj_filt, "otu_table") # call this by the otu-table
taxa_obj_filt$data$tax_proportions <- calc_obs_props(taxa_obj_filt, "tax_data")
########################################################################################

# i don't yet understand why the calc_taxon_abund() step is dropping us down to 755, when all of our taxa are unique already... but oh well
pseq_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_obj_filt, "tax_proportions", cols = target_samples_no_neg)

taxa_obj_filt$data$tax_abund <- calc_taxon_abund(taxa_obj_filt, "tax_proportions", cols = target_samples_no_neg)

# not sure if we want the groups argument here or not
pseq_obj_filt$data$tax_occ <- calc_n_samples(pseq_obj_filt, "tax_abund")
taxa_obj_filt$data$tax_occ <- calc_n_samples(taxa_obj_filt, "tax_abund")

# doing the compare now with the abund table
pseq_obj_filt$data$compare_tax_abund <- compare_groups(pseq_obj_filt,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_obj_filt$data$sample_data$case)



taxa_obj_filt$data$compare_tax_abund <- compare_groups(taxa_obj_filt,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups =pseq_obj_filt$data$sample_data$case)



range(pseq_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)

range(taxa_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)

taxa_obj_filt_slim<-taxa_obj_filt$data$compare_tax_abund
taxa_obj_filt$data$compare_tax_abund
taxa_obj_filt$data$compare_tax_abund%>%filter(taxon_id=="bow")

pseq_obj_filt$data$tax_data%>%filter(genus=="Sphingomonas")

pseq_obj_filt$data$compare_tax_abund%>%filter(taxon_id=="aoz")
#%>%filter(wilcox_p_value<0.05) 
df2<-taxa_obj_filt$data$compare_tax_abund%>%filter(treatment_1=="COVID19"|treatment_2=="Uninfected")%>%filter(wilcox_p_value<=0.08)#%>%filter(abs(log2_median_ratio)>=1.76)
df3<-taxa_obj_filt$data$compare_tax_abund%>%filter(treatment_1=="COVID19")%>%filter(treatment_2=="CAP")%>%filter(wilcox_p_value<=0.9)#%>%filter(abs(log2_median_ratio)>=1.76)
#df31<-taxa_obj_filt$data$compare_tax_abund%>%filter(treatment_1=="Uninfected")%>%filter(treatment_2=="CAP")%>%filter(wilcox_p_value<=0.3)#%>%filter(abs(log2_median_ratio)>=1.76)

df2
df3
df31
df31
df4<-full_join(df2,df3)
#df4<-full_join(df31,df4)
df4


df5<-inner_join(df4,pseq_obj_filt$data$tax_data)%>%filter(!is.na(genus))%>%filter(wilcox_p_value<=0.05)
df5





df6<-df5%>%filter(!is.na(species))%>%distinct_all()%>%filter(duplicated(genus))%>%arrange(treatment_1,treatment_2)
df6
as.data.frame(df6)
df4
hist(df4$wilcox_p_value)
keep<-df4$taxon_id
pseq_obj
pseq_obj_filt<-pseq_obj%>%filter_taxa(taxon_ids%in%keep,supertaxa = T)
pseq_obj_filt
taxa_obj_filt

##########re-run the tax stuff#########
pseq_obj_filt<-pseq_obj%>%filter_taxa(taxon_ids%in%keep,supertaxa = T)
pseq_obj_filt
pseq_obj_filt$data$tax_proportions <- calc_obs_props(pseq_obj_filt, "otu_table") # call this by the otu-table
pseq_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_obj_filt, "tax_proportions", cols = target_samples_no_neg)
pseq_obj_filt$data$tax_occ <- calc_n_samples(pseq_obj_filt, "tax_abund")
pseq_obj_filt$data$compare_tax_abund <- compare_groups(pseq_obj_filt,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_obj_filt$data$sample_data$case)
warnings()
range(pseq_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
pseq_obj_filt$data$compare_tax_abund$log2_median_ratio
taxon_ids(obj = pseq_obj_filt)
pseq_obj_filt
# we're looking at -6.8 to 7.3

# trying a plot (commenting out node_label = taxon_names makes them plot much faster once saved to an object, good for when trying different things)
# 
# fig1 <- heat_tree_matrix(taxa_obj_filt,
#                          data = "compare_tax_abund",
#                          node_size = n_obs,
#                          node_label = taxon_names,
#                          node_color = log2_median_ratio,
#                          node_color_range = diverging_palette(),
#                          node_color_trans = "linear",
#                          node_color_interval = c(-3,3),
#                          node_size_axis_label = "Taxa counts",
#                          node_color_axis_label = "Log2 median ratio",
#                          layout = "davidson-harel",
#                          initial_layout = "reingold-tilford",
#                          output_file = "node_color_trans_linear.pdf",
#                          verbose = TRUE)
# fig1
range(pseq_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
pseq_fig1 <- heat_tree_matrix(pseq_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-2,2),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         key_size = 0.7,
                         seed = 314,  
                         aspect_ratio = 1.5,
                         initial_layout = "reingold-tilford",
                         repel_labels=T,
                         output_file = "node_color_trans_linear.pdf",
                         verbose = TRUE)
pseq_fig1

# 
# 
# fig2 <- heat_tree_matrix(taxa_obj_filt,
#                          data = "compare_tax_abund",
#                          node_size = n_obs,
#                          node_label = taxon_names,
#                          node_color = log2_median_ratio,
#                          node_color_range = diverging_palette(),
#                          node_color_trans = "area",
#                          node_color_interval = c(-9,9),
# #                         node_color_interval = c(-7.5,7.5),
#                          node_size_axis_label = "Taxa counts",
#                          node_color_axis_label = "Log2 median ratio",
#                          layout = "davidson-harel",
#                          initial_layout = "reingold-tilford",
#                          output_file = "node_color_trans_area.pdf",
#                          key_size = 0.5,
#                          seed = 314,  
#                          title = "Taxonomic Comparison of COVID19 vs Uninfected & CAP",
#                          title_size = 0.08,
#                          verbose = TRUE)
# fig2
pseq_fig2 <- heat_tree_matrix(pseq_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-3.9,3.9),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "node_color_trans_area.pdf",
                         verbose = TRUE)





#alpha diversity

library(vegan)

pseq_obj_filt$data$sample_data$inv_simp <- diversity(pseq_obj_filt$data$otu_table[, pseq_obj_filt$data$sample_data$sample_id],
                                                     index = "invsimpson",
                                                     MARGIN = 2) # What orietation the matrix is in
library(ggplot2)
ggplot(pseq_obj_filt$data$sample_data, aes(x = case, y = inv_simp)) +
  geom_boxplot()
anova_result <- aov(inv_simp ~ case, pseq_obj_filt$data$sample_data)
summary(anova_result)
# Df Sum Sq Mean Sq F value Pr(>F)  
# case         2    992   495.8    2.83 0.0647 .
# Residuals   83  14537   175.1                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Tukey’s Honest Significant Difference (HSD) 
library(agricolae)
tukey_result <- HSD.test(anova_result, "case", group = TRUE)
print(tukey_result)
# $statistics
# MSerror Df     Mean       CV
# 175.1488 83 15.48333 85.47501
# 
# $parameters
# test name.t ntr StudentizedRange alpha
# Tukey   case   3         3.374985  0.05
# 
# $means
# inv_simp      std  r      Min      Max      Q25       Q50      Q75
# Control_Healthy 20.04256 17.78298 29 2.580625 70.27482 8.188479 14.876679 23.64051
# Control_Sick    14.52070 10.11376 25 1.218328 42.50921 6.817982 14.506570 17.72547
# COVID19         12.10359 10.20416 32 1.791336 46.79076 7.475984  9.352815 13.03280
# 
# $comparison
# NULL
# 
# $groups
# inv_simp groups
# Control_Healthy 20.04256      a
# Control_Sick    14.52070      a
# COVID19         12.10359      a
# 
# attr(,"class")
# [1] "group"


group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggplot(pseq_obj_filt$data$sample_data, aes(x = case, y = inv_simp)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data), y = max(pseq_obj_filt$data$sample_data$inv_simp) + 3, label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot() +
  ggtitle("Alpha diversity of disease case sites") +
  xlab("disease case") +
  ylab("Inverse Simpson Index")
library(ggpubr)
ggboxplot(data =pseq_obj_filt$data$sample_data,x = "case",y = "inv_simp",add = "jitter",color="case")+
  stat_compare_means(method = "t.test", ref.group = "Control_Healthy")

+stat_compare_means()
