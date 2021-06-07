library(tidyverse)
library(phyloseq)
setwd("c:/github/microbial/step2_kraken2_analysis")


# if wanting to load the objects
#load("heat-tree.RData")

# reading in metadata
sample_info_tab <- read.table("Combined_BALF_GO_Terms_metadata.tsv", header = TRUE, sep = "\t")
sample_info_tab$outcome<-gsub("stabilized","Survived",x =sample_info_tab$outcome)
sample_info_tab$outcome<-gsub("recovered","Survived",x =sample_info_tab$outcome)
sample_info_tab$outcome<-gsub("Stabilized","Survived",x =sample_info_tab$outcome)
sample_info_tab$outcome<-gsub("Recovered","Survived",x =sample_info_tab$outcome)
sample_info_tab$outcome<-gsub("deceased","Deceased",x =sample_info_tab$outcome)
unique(sample_info_tab$outcome)
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
write.table(perc_remain,"percent_remaining_post_decontam.tsv",sep = "\t",quote = F,row.names = F)
tax<-as_tibble(bac_tax_info_tab,rownames="taxid")
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


pseq_decontam<-merge_phyloseq(pseq_not_shen,pseq_shen_decontam)
pseq_decontam_no_neg<-subset_samples(physeq = pseq_decontam,case!="Control_Neg")

pseq_decontam_no_neg_core<-core(x = pseq_decontam_no_neg,detection = 0,prevalence = 50/100)
pseq_decontam_no_neg_core #[ 566 taxa and 86 samples ]



outcome_pseq_decontam_no_neg_core<-subset_samples(pseq_decontam_no_neg_core,case=="COVID19")
outcome_pseq_decontam_no_neg_core<-subset_samples(outcome_pseq_decontam_no_neg_core,outcome!="NA")
range(taxa_sums(outcome_pseq_decontam_no_neg_core))



# 
 library(metacoder)
# # making our data look like theirs in this example page: https://grunwaldlab.github.io/metacoder_documentation/example.html
#   # which we can look at in the hmp_otus object after loading the metacoder library
# hmp_otus
#    # they have an ID column, a lineage column, then the counts, we're going to turn our tax info into their lineage column format
# hmp_otus$lineage %>% head(30) %>% tail
# 
# lineage_vector <- paste0(
#     "d__", bac_tax_info_tab$domain, ";",
#     "p__", bac_tax_info_tab$phylum, ";",
#     "c__", bac_tax_info_tab$class, ";",
#     "o__", bac_tax_info_tab$order, ";",
#     "f__", bac_tax_info_tab$family, ";",
#     "g__", bac_tax_info_tab$genus, ";",
#     "s__", bac_tax_info_tab$species)
# 
# # and the don't have the "NA" ones, they just have nothing for them, so removing that business so ours stop at the last rank present like theirs
# lineage_vector <- gsub("(;[c-s]__NA)*$", "", lineage_vector)
# head(lineage_vector)
# # beautiful, building table like hmp_otus
# 
# # outcome_tax_counts_tab<-data.frame(otu_table(outcome_pseq_decontam_no_neg_core))
# # outcome_tax_info_tab<-data.frame(tax_table(outcome_pseq_decontam_no_neg_core))
# tab_for_taxa_object <- data.frame(taxid = row.names(bac_tax_info_tab), lineage=lineage_vector, bac_tax_counts_tab)

# #lets get rid of the non covid19 outcome samples
# keep_vector<-rownames(meta(outcome_pseq_decontam_no_neg_core))
# keep_vector
# tab_for_taxa_object<-tab_for_taxa_object%>%pivot_longer(-c(lineage,taxid),names_to = "SRA")%>%filter(SRA %in% keep_vector)%>%pivot_wider(names_from=SRA,values_from = value)
# head(tab_for_taxa_object)
# head(hmp_otus)[, 1:10]
# str(tab_for_taxa_object)
# i think all's ok, let's give it a shot

pseq_obj<-parse_phyloseq(outcome_pseq_decontam_no_neg_core)
print(pseq_obj)

pseq_obj_filt <- taxa::filter_taxa(pseq_obj, n_obs > 5)
pseq_obj_filt$data$tax_data <- zero_low_counts(pseq_obj_filt, data = "tax_data", min_count = 5)



no_reads <- rowSums(pseq_obj_filt$data$otu_table[, pseq_obj_filt$data$sample_data$sample_id]) == 0

sum(no_reads)
pseq_obj_filt <- filter_obs(pseq_obj_filt, data = "otu_table", ! no_reads, drop_taxa = TRUE)
pseq_obj_filt

# this is where we normalize across samples (i'm curious what this spits out when there are negative values in the table like from the vst, ha)
#################################################################################
pseq_obj_filt$data$tax_proportions <- calc_obs_props(pseq_obj_filt, "otu_table")

########################################################################################

# i don't yet understand why the calc_taxon_abund() step is dropping us down to 755, when all of our taxa are unique already... but oh well

pseq_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_obj_filt, "tax_proportions", cols = pseq_obj_filt$data$sample_data$sample_id)

# not sure if we want the groups argument here or not
pseq_obj_filt$data$tax_occ <- calc_n_samples(pseq_obj_filt, 
                                             "tax_abund", 
                                             cols = pseq_obj_filt$data$sample_data$sample_id,
                                             groups = pseq_obj_filt$data$sample_data$outcome)



pseq_obj_filt$data$tax_occ 

set.seed(314) # This makes the plot appear the same each time it is run 
heat_tree(pseq_obj_filt, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Deceased, 
          node_size_axis_label = "Read Counts",
          node_color_axis_label = "Samples with reads",
          title="Deceased",
          repel_labels=T,
          aspect_ratio=0.95,
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations
heat_tree(pseq_obj_filt, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Survived, 
          repel_labels=T,
          aspect_ratio=0.95,
          node_size_axis_label = "Read Counts",
          node_color_axis_label = "Samples with reads",
          title="Survived",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford")


pseq_obj_filt$data$diff_table <- compare_groups(pseq_obj_filt,
                                      dataset = "tax_abund",
                                      cols = pseq_obj_filt$data$sample_data$sample_id,
                                      groups = pseq_obj_filt$data$sample_data$outcome,
                                      combinations = list(c('Deceased', 'Survived')))

range(pseq_obj_filt$data$diff_table$log2_median_ratio)

set.seed(123)
heat_tree(pseq_obj_filt, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-5, 5), # The range of `log2_median_ratio` to display
          node_color_range = c("#00A1D5FF","grey75","#B24745FF"), # The color palette used
          node_size_axis_label = "OTU count",
          # node_size_range = c(0.01, 0.1),
          # edge_size_range = c(0.005, 0.01),
          node_color_axis_label = "Log 2 ratio of median proportions",
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          verbose=T,
          title="Deceased vs Survived",
          title_size=0.05) # The layout algorithm that initializes node locations




pseq_obj_filt$data$diff_table$p.adjust<-p.adjust(p = pseq_obj_filt$data$diff_table$wilcox_p_value, method = "BH", n = 84)
pseq_obj_filt$data$diff_table$p.adjust

pseq_obj_filt$data$diff_table_sig<-pseq_obj_filt$data$diff_table%>%filter(wilcox_p_value<0.05)

pseq_sig<-as_tibble(pseq_obj_filt$data$diff_table_sig)
names<-as.data.frame(pseq_obj_filt$taxon_names())
names$taxon_id<-rownames(names)
names
sig<-inner_join(pseq_sig,names)
sig





#alpha diversity

library(vegan)

pseq_obj_filt$data$sample_data$inv_simp <- diversity(pseq_obj_filt$data$otu_table[, pseq_obj_filt$data$sample_data$sample_id],
                                  index = "invsimpson",
                                  MARGIN = 2) # What orietation the matrix is in
library(ggplot2)
ggplot(pseq_obj_filt$data$sample_data, aes(x = outcome, y = inv_simp)) +
  geom_boxplot()
anova_result <- aov(inv_simp ~ outcome, pseq_obj_filt$data$sample_data)
summary(anova_result)
# Df Sum Sq Mean Sq F value Pr(>F)
# outcome      1  124.5   124.5   1.041  0.318
# Residuals   23 2750.6   119.6


#Tukey’s Honest Significant Difference (HSD) 
library(agricolae)
tukey_result <- HSD.test(anova_result, "outcome", group = TRUE)
print(tukey_result)
# 
# $statistics
# MSerror Df     Mean       CV
# 119.591 23 13.52398 80.86204
# 
# $parameters
# test  name.t ntr StudentizedRange alpha
# Tukey outcome   2         2.925523  0.05
# 
# $means
# inv_simp       std  r      Min      Max      Q25       Q50      Q75
# Deceased 10.79122  2.661937 10 7.939966 14.55147 8.692691  9.529152 13.17032
# Survived 15.34582 13.853363 15 2.272643 46.79076 7.671496 10.953091 12.96878
# 
# $comparison
# NULL
# 
# $groups
# inv_simp groups
# Survived 15.34582      a
# Deceased 10.79122      a
# 
# attr(,"class")
# [1] "group"



group_data <- tukey_result$groups[order(rownames(tukey_result$groups)),]
ggplot(pseq_obj_filt$data$sample_data, aes(x = outcome, y = inv_simp)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data), y = max(pseq_obj_filt$data$sample_data$inv_simp) + 3, label = group_data$groups),
            col = 'black',
            size = 10) +
  geom_boxplot() +
  ggtitle("Alpha diversity of disease outcome sites") +
  xlab("disease outcome") +
  ylab("Inverse Simpson Index")
