################################################################################
################################################################################
# COVID19 Disease Outcome Metacoder Heat Tree Data Visualiz(s)ations
# Jochum, Michael D.
# 19 May 2021
################################################################################
################################################################################

vst_physeq 


library(tidyverse)
library(phyloseq)
setwd("~/github/microbial/step2_kraken2_analysis")


# if wanting to load the objects
#load("outcome_heat_tree_data.RData")


################################################################################
#preprocessing metadata, taxonomy and read count files
################################################################################
# reading in metadata
sample_info_tab <- read.table("Combined_BALF_GO_Terms_metadata.tsv", header = TRUE, sep = "\t")
# reading in tax file
full_tax_tab <- read.table("all_step2_kraken2_summaries.tsv", header = TRUE, row.names = 1, sep = "\t", quote = '', comment.char="")

# making table of just counts (taxid as rows)
tax_counts_tab <- full_tax_tab[, grep("counts", names(full_tax_tab))]
# shortening names
names(tax_counts_tab) <- sub("_step2_kraken2_mason_db_read_counts", "", names(tax_counts_tab))
# # making table of just tax info
tax_info_tab <- full_tax_tab[, 1:7]
# # filtering down to only the 86 samples in the GO term analysis
target_samples <- scan("GO_TERM_samples.tsv", what = "character")[-1]
target_samples_no_neg<-target_samples
# #add the negative controls back into the analysis
negs<-c("CRR125995",'CRR125996',"CRR125997","CRR125998")
target_samples<-c(target_samples,negs)
#length(target_samples) # 90
#
target_tax_counts_tab <- tax_counts_tab %>% select(all_of(target_samples))
#dim(target_tax_counts_tab) #11044    90

# another look they were all found even though right size table came through
#setdiff(target_samples, names(target_tax_counts_tab))
#
target_sample_info_tab <- sample_info_tab %>% filter(accession %in% target_samples)
#dim(target_sample_info_tab)#[1] 90 74
#setdiff(target_samples, target_sample_info_tab$accession)
#
target_sample_info_tab <- target_sample_info_tab %>% column_to_rownames("accession")
# filtering down to just Bacteria
bac_tax_counts_tab <- target_tax_counts_tab[!is.na(tax_info_tab$domain) & tax_info_tab$domain == "Bacteria", ]
# filter out only bacterial reads for the figure
bac_tax_info_tab <- tax_info_tab[!is.na(tax_info_tab$domain) & tax_info_tab$domain == "Bacteria", ] #[1] 8661    7
dim(bac_tax_info_tab)
#all.equal(row.names(bac_tax_info_tab), row.names(bac_tax_counts_tab))


################################################################################
#make a phyloseq object
################################################################################
#MJ-Make a pseq-object so you can have something to filter in and out of for decontam and metacoder)
pseq<-phyloseq(otu_table(bac_tax_counts_tab,taxa_are_rows = T),
               tax_table(as.matrix(bac_tax_info_tab)),
               sample_data(target_sample_info_tab))
pseq

################################################################################
#decontaminate the reads with negative conrols
################################################################################
library(decontam)
#make a pseq obj with just the shen samples that have negative controls
pseq_shen<-subset_samples(pseq,bioproject=="CRA002476")
pseq_shen #[ 8661 taxa and 49 samples ]
# make a pseq obj with all the other samples for merging later on
pseq_not_shen<-subset_samples(pseq,bioproject!="CRA002476") #[ 8661 taxa and 49 samples ]
pseq_not_shen #  [ 8661 taxa and 41 samples ]


library(microbiome)
# export the tax table to a data frame
tax_df_shen<-data.frame(tax_table(pseq_shen))
tax_df_shen
#export the metadata to a dataframe a make a vector of the negative controls
sam_shen<-meta(pseq_shen)
vector_for_decontam <- sam_shen$case=="Control_Neg"
vector_for_decontam                    


count_df_shen<-as.data.frame(otu_table(pseq_shen))
colSums(count_df_shen)

#itdenitfy the likely contaminants 
contam_df_shen <- isContaminant(t(count_df_shen), neg=vector_for_decontam)
table(contam_df_shen$contaminant) 
# FALSE  TRUE 
#4134  4527 

# getting vector holding the identified contaminant IDs
contam_shen <- row.names(contam_df_shen[contam_df_shen$contaminant == TRUE, ])
contam_shen<-as.factor(contam_shen)
contam_shen



count_df_shen<-as_tibble(count_df_shen,rownames = "taxid")
count_df_shen
#count_df_shen[1]<-NULL
#filter out the contaminants
clean_count_df_shen<-count_df_shen%>%filter(!taxid%in%contam_shen)

a<-clean_count_df_shen%>%select(-taxid)
b<-count_df_shen%>%select(-taxid)
perc_remain<-(colSums(a)/colSums(b))*100
mean(perc_remain) # 88.86533 we kept about 88% of our reads (I think)
write.table(perc_remain,"percent_remaining_post_decontam.tsv",sep = "\t",quote = F,row.names = F)
tax<-as_tibble(bac_tax_info_tab,rownames="taxid")
#make tables of the contaminants and not conaminants tax_tables
contaminants<-tax%>%filter(taxid %in% contam_shen)
not_contaminants<-tax%>%filter(!taxid %in% contam_shen)
#write out a table of each of the outputs
write.table(x = contaminants,file = "Shen_contaminants.tsv",sep = "\t",row.names = F,quote = F)
write.table(x = not_contaminants,file = "Shen_not_contaminants.tsv",sep = "\t",row.names = F,quote = F)

#pull out the contaminants from the shen pseq obj
pseq_shen_decontam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)
pseq_shen_decontam

#merge the decontaminated shen phyloseq objects with the rest of the samples
pseq_decontam<-merge_phyloseq(pseq_not_shen,pseq_shen_decontam)
#remove the negative controls from the analysis
pseq_decontam_no_neg<-subset_samples(physeq = pseq_decontam,case!="Control_Neg")

#make a 50% prevalence and detection cutoff
pseq_decontam_no_neg_core<-core(x = pseq_decontam_no_neg,detection = 50/100,prevalence = 50/100)
pseq_decontam_no_neg_core #[ 566 taxa and 86 samples ]


#subset the phyloseq object with just the COVID19 samples with known survival outcomes
survival_pseq_decontam_no_neg_core<-subset_samples(pseq_decontam_no_neg_core,case=="COVID19")
survival_pseq_decontam_no_neg_core<-subset_samples(pseq_decontam_no_neg_core,survival!="NA")

survival_pseq_decontam_no_neg_core #[ 566 taxa and 25 samples ]
range(taxa_sums(survival_pseq_decontam_no_neg_core)) #124 3071073
range(sample_sums(survival_pseq_decontam_no_neg_core)) # 267 9385857


################################################################################
############## Metacoder hat tree data visualiz(s)ations #######################
################################################################################
library(metacoder)

pseq_obj<-parse_phyloseq(survival_pseq_decontam_no_neg_core)
pseq_obj

#convert all counts below a minimum number (5) to zero
pseq_obj$data$tax_data <- zero_low_counts(pseq_obj, data = "tax_data", min_count = 5)
#filter taxa with less than 5 observations
pseq_obj_filt <- taxa::filter_taxa(pseq_obj, n_obs > 5)
pseq_obj_filt

# #double check to make sure there arent any rows with zero counts remaining in the dataset
 no_reads <- rowSums(pseq_obj_filt$data$otu_table[, pseq_obj_filt$data$sample_data$sample_id]) == 0
 sum(no_reads) # [1] 0 good
# 
# #just in case there were, we could remove them by filtering the observations
 pseq_obj_filt <- filter_obs(pseq_obj_filt, data = "otu_table", ! no_reads, drop_taxa = TRUE)
 pseq_obj_filt

# Normalize observations across samples 
#################################################################################
pseq_obj_filt$data$tax_proportions <- calc_obs_props(pseq_obj_filt, "otu_table")
#Calculating proportions from counts for 25 columns for 566 observations.
#figure out the range in the normalized read counts
#ps dont use the character taxid column 1
range(pseq_obj_filt$data$tax_proportions[-1])#[1] 0.0000000 0.6519896 
################################################################################

#Sum up the observation values for each taxon
pseq_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_obj_filt, 
                                                 "tax_proportions", 
                                                 cols = pseq_obj_filt$data$sample_data$sample_id)
#Summing per-taxon counts from 25 columns for 84 taxa
range(pseq_obj_filt$data$tax_abund[-1]) # 0 1 ... good

# Count the number of samples
pseq_obj_filt$data$tax_occ <- calc_n_samples(pseq_obj_filt, 
                                             "tax_abund", 
                                             cols = pseq_obj_filt$data$sample_data$sample_id,
                                             groups = pseq_obj_filt$data$sample_data$survival)

#Calculating number of samples with a value greater than 0 for 25 columns in 2 groups for 84 observations
pseq_obj_filt$data$tax_occ 

#make the individual data visualizations for each treatment
set.seed(314) # This makes the plot appear the same each time it is run 

# deceased heat tree
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

#survived heat tree
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

################################################################################
#Compare the abundance of the taxon between deceased and survived
################################################################################
pseq_obj_filt$data$diff_table <- compare_groups(pseq_obj_filt,
                                      data = "tax_abund",
                                      cols = pseq_obj_filt$data$sample_data$sample_id,
                                      groups = pseq_obj_filt$data$sample_data$survival,
                                      combinations =list(c('Survived', 'Deceased')))

range(pseq_obj_filt$data$diff_table$log2_median_ratio) #[1] -5.175576  5.207248
################################################################################

#make the visualization
set.seed(314)
heat_tree(pseq_obj_filt, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-5.175576 , 5.207248), # The range of `log2_median_ratio` to display
          node_color_range = c("#00A1D5FF","grey75","#B24745FF"), # The color palette used
          node_size_axis_label = "Read count",
          # node_size_range = c(0.01, 0.1),
          # edge_size_range = c(0.005, 0.01),
          node_color_axis_label = "Log 2 ratio of median proportions",
          initial_layout = "reingold-tilford",
          layout = "davidson-harel",
          verbose=T,
          title="Deceased vs Survived",
          title_size=0.05) # The layout algorithm that initializes node locations

#adjust the p value for multiple test correction using n = 84 taxon comparisons 
pseq_obj_filt$data$diff_table$p.adjust<-p.adjust(p = pseq_obj_filt$data$diff_table$wilcox_p_value, method = "BH", n = 84)
pseq_obj_filt$data$diff_table$p.adjust

#make a table with the statistically significant taxon ids 
pseq_obj_filt$data$diff_table_sig<-as_tibble(pseq_obj_filt$data$diff_table%>%filter(wilcox_p_value<0.05))
#get a list of the taxonomy names that match the taxon ids
names<-as.data.frame(pseq_obj_filt$taxon_names())
names$taxon_id<-rownames(names)
#merge the two datasets together
sig<-inner_join(pseq_obj_filt$data$diff_table_sig,names)
sig
# # A tibble: 14 x 9
# taxon_id treatment_1 treatment_2 log2_median_ratio median_diff mean_diff wilcox_p_value p.adjust `pseq_obj_filt$taxon_names()`
# <chr>    <chr>       <chr>                   <dbl>       <dbl>     <dbl>          <dbl>    <dbl> <chr>                        
#   1 aat      Deceased    Survived                 2.29     0.425     0.296         0.0264    0.185   Betaproteobacteria           
# 2 aaz      Deceased    Survived                -5.13    -0.103    -0.104         0.0308    0.199   Bacteroidia                  
# 3 abd      Deceased    Survived                 1.84     0.0549    0.130         0.0137    0.124   Bacilli                      
# 4 acn      Deceased    Survived                 2.24     0.403     0.297         0.0163    0.124   Burkholderiales              
# 5 act      Deceased    Survived                 2.97     0.00242   0.00210       0.00353   0.0740  Vibrionales                  
# 6 acw      Deceased    Survived                -5.18    -0.0998   -0.102         0.00962   0.124   Bacteroidales                
# 7 ads      Deceased    Survived                 3.16     0.00183   0.00171       0.0157    0.124   Alteromonadales              
# 8 agd      Deceased    Survived                 2.25     0.361     0.371         0.000165  0.00691 Comamonadaceae               
# 9 ags      Deceased    Survived                 2.97     0.00242   0.00210       0.00353   0.0740  Vibrionaceae                 
# 10 ahe      Deceased    Survived                 1.77     0.0111    0.0645        0.0475    0.274   Streptococcaceae             
# 11 aii      Deceased    Survived                 3.61     0.00361   0.00387       0.0156    0.124   Yersiniaceae                 
# 12 anj      Deceased    Survived                 2.10     0.00470   0.00435       0.0156    0.124   Salmonella                   
# 13 ann      Deceased    Survived                 3.80     0.00169   0.00181       0.00492   0.0827  Vibrio                       
# 14 apt      Deceased    Survived                 5.21     0.405     0.377         0.000165  0.00691 Variovorax   

#save.image(file = "outcome_heat_tree_data.RData")

################################################################################
#Exit 00
