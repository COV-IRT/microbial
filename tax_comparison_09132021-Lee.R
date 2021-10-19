library(tidyverse)
library(phyloseq)
library(Maaslin2)
################################################################################
##################### Data import / Preprocessing ##############################
################################################################################
getwd()
#setwd("c:/github/microbial/step2_kraken2_analysis/")

# if wanting to load the objects
#load("heat-tree-Jochum.RData")

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


################################################################################
##################### decontam ##############################
################################################################################

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

#quicky regexy stuff to update the case terms
sample_data(pseq_decontam_no_neg)$case<-gsub('Control_Healthy','Uninfected',sample_data(pseq_decontam_no_neg)$case)
sample_data(pseq_decontam_no_neg)$case<-gsub('COVID 19','COVID-19',sample_data(pseq_decontam_no_neg)$case)
sample_data(pseq_decontam_no_neg)$case<-gsub('Control_Sick','CAP',sample_data(pseq_decontam_no_neg)$case)

################################################################################
################################ Maaslin2 ######################################
################################################################################

library(Maaslin2)

#make a data frame containine the abundances and metadata for the counts
df_input_data<-data.frame(t(otu_table(pseq_decontam_no_neg_core)))
df_input_metadata<-data.frame(sample_data(pseq_decontam_no_neg_core))

#run Maaslin2 on the taxa counts using following paramenters

# CLR normalization
# Log transformation
# 0.1 and 0.1 abundance and prevalence cutoffs
# Bejamini Hochberg multiple test correction

#########################################
# warning this takes FOREVER to run
##########################################
case_norm<-Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output="./terms_vs_case_comp_norm_09102021",
  min_abundance = 0.01,
  min_prevalence = 0.1,
  normalization = "CLR",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("sample_name","publication"),
  fixed_effects = c("case"),
  correction="BH",
  standardize = TRUE,
  cores = 20,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n =100,
  reference="case,COVID19")

###########################################

#make a copy of the results
res<-as_tibble(case_norm$results)
#filter the results to only include significant taxa
keep<-res%>%filter(pval<0.05)
#fix the character
keep<-keep%>%mutate(feature=gsub("X","",feature))

#make a phyloseq_object that has subsets only the stat. significant taxa
physeq_keep<-prune_taxa(taxa = keep$feature, x =  pseq_decontam_no_neg)

physeq_keep #[ 516 taxa and 86 samples ]

# make a key containing the taxonomy of the sig. taxa
keep_tax<-data.frame(tax_table(physeq_keep))
keep_tax$feature<-rownames(keep_tax)

keep<-inner_join(keep,keep_tax)
Sphing<-keep%>%filter(genus=="Sphingomonas")
Sphing

length(unique(keep$genus))

# make a hack that recreates the metaphlan style lineage
keep<-keep%>%mutate(lineage=paste0("d__",keep$domain,
                                   ";p__",keep$phylum,
                                   ";c__",keep$class,
                                   ";o__",keep$order,
                                   ";f__",keep$family,
                                   ";g__",keep$genus,
                                   ";s__",keep$species),
                    lineage=gsub("d__NA;p__NA;c__NA;o__NA;f__NA;g__NA;s__NA",NA,lineage),
                    lineage=gsub(";p__NA;c__NA;o__NA;f__NA;g__NA;s__NA","",lineage),
                    lineage=gsub(";c__NA;o__NA;f__NA;g__NA;s__NA","",lineage),
                    lineage=gsub(";o__NA;f__NA;g__NA;s__NA","",lineage),
                    lineage=gsub(";f__NA;g__NA;s__NA","",lineage),
                    lineage=gsub(";g__NA;s__NA","",lineage),
                    lineage=gsub(";s__NA","",lineage))

head(keep$lineage) # looks good

############# MIKE LEE ADDITIONS 14-Sep-2021 ##############
#load("heat-tree-Jochum.RData")
#filter the results to only include significant taxa
keep <- res %>% filter(qval < 0.05)

#fix the character
keep <- keep %>% mutate(feature = gsub("X", "", feature))
keep
vector_of_OTU_IDs_we_want <- keep %>% pull(feature)


length(vector_of_OTU_IDs_we_want)
keep
############################################################

############# MIKE LEE ADDITIONS 14-Sep-2021 ##############
   # making what I think is a full one following what was done above to the filtered ones, but the process of subsetting and replacing the values will be the same when you do it with what you know to be the wanted full one
library(tidyverse)

full_pseq_obj <- parse_phyloseq(pseq_decontam_no_neg)
full_pseq_obj$data$tax_proportions <- calc_obs_props(full_pseq_obj, "otu_table")
full_pseq_obj$data$tax_abund <- calc_taxon_abund(full_pseq_obj, "tax_proportions", cols = target_samples_no_neg)
full_pseq_obj$data$tax_occ <- calc_n_samples(full_pseq_obj, "tax_abund")
full_pseq_obj$data$compare_tax_abund <- compare_groups(full_pseq_obj,
                                                  "tax_abund",
                                                  cols = target_samples_no_neg,
                                                  groups = full_pseq_obj$data$sample_data$case)

full_pseq_obj



# the full compare table has taxon_ids
full_pseq_obj$data$compare_tax_abund
# we have OTU IDs of what we want
vector_of_OTU_IDs_we_want
# these are mapped in this table
full_pseq_obj$data$tax_data
# so getting the taxon_ids we want from that table based on the OTU IDs we have
vector_of_wanted_taxon_ids <- full_pseq_obj$data$tax_data %>% filter(otu_id %in% vector_of_OTU_IDs_we_want) %>% pull(taxon_id) %>% unique()
length(vector_of_OTU_IDs_we_want) # 161
length(vector_of_wanted_taxon_ids) # 128 # this is going to be shorter because some OTUs can belong to the same taxon_id

# now making subset compare_tax_abund table holding just the taxon_ids we want
wanted_compare_tab <- full_pseq_obj$data$compare_tax_abund %>% filter(taxon_id %in% vector_of_wanted_taxon_ids)
dim(wanted_compare_tab) # 384 # why the f is this larger than 128?
head(wanted_compare_tab)
data.frame(table(wanted_compare_tab %>% pull(taxon_id)))

wanted_compare_tab %>% filter(taxon_id == "aab") # oh yeah, it has to have values for each contrast, which is why it's 384 (128 * 3)

# making copy of object we are modifying
mod_pseq_obj <- full_pseq_obj

# filtering that down to the size we want
filt_mod_pseq_obj <- mod_pseq_obj %>% filter_taxa(taxon_ids %in% vector_of_wanted_taxon_ids, reassign_obs = FALSE, drop_obs = TRUE) # gave a warning about "dataset 3", but i can't figure out what that means and it doesn't seem to affect anything we want

# replacing compare_tax_abund table (before doing so, it's the same as before we filtered)
filt_mod_pseq_obj$data$compare_tax_abund <- as_tibble(wanted_compare_tab)

# now this is in here with the same info from the full compare_groups run
filt_mod_pseq_obj$data$compare_tax_abund
full_pseq_obj$data$compare_tax_abund
# here's looking at taxon_id "aci"
filt_mod_pseq_obj$data$compare_tax_abund %>% filter(taxon_id == "aci")
full_pseq_obj$data$compare_tax_abund %>% filter(taxon_id == "aci")

# # making heattree
# filt_mod_pseq_fig1 <- heat_tree_matrix(filt_mod_pseq_obj,
#                               data = "compare_tax_abund",
#                               node_size = n_obs,
#                               node_label = taxon_names,
#                               node_color = log2_median_ratio,
#                               node_color_range = diverging_palette(),
#                               node_color_trans = "area",
#                               node_color_interval = range(filt_mod_pseq_obj$data$compare_tax_abund$log2_median_ratio, finite = TRUE),
#                               node_size_axis_label = "Taxa counts",
#                               node_color_axis_label = "Log2 median ratio",
#                               layout = "davidson-harel",
#                               key_size = 0.6,
#                               seed = 314,
#                               aspect_ratio = 1,
#                               initial_layout = "reingold-tilford",
#                               repel_labels=T,
#                               output_file = "pseq_obj_filt_node_color_trans_linear.pdf",
#                               verbose = TRUE)

# filt_mod_pseq_fig1

## this looks a little odd, I think because we eliminated all the intermediate taxa nodes
## like if we got rid of higher ranks but kept a species rank if it was significant, so the
## branch just juts straight out from the domain
vector_of_wanted_taxon_ids

### going to try doing it a little differently (filtering keeping the parent nodes, then adding those to what we pull as our wanted table)
# filtering that down to the size we want
filt_mod_pseq_obj_with_parents <- mod_pseq_obj %>% filter_taxa(taxon_ids %in% vector_of_wanted_taxon_ids, supertaxa = TRUE, reassign_obs = FALSE, drop_obs = TRUE)

# replacing compare_tax_abund table (before doing so, it's the same as before)
    # this time, we need to remake the wanted table, including the parent taxon_ids
vector_of_wanted_taxon_ids_with_parents <- filt_mod_pseq_obj_with_parents$taxa %>% names
length(vector_of_wanted_taxon_ids_with_parents)

wanted_compare_tab_with_parents <- mod_pseq_obj$data$compare_tax_abund %>% filter(taxon_id %in% vector_of_wanted_taxon_ids_with_parents)

## replacing compare tab
filt_mod_pseq_obj_with_parents$data$compare_tax_abund <- as_tibble(wanted_compare_tab_with_parents)


# making heattree
filt_mod_pseq_with_parents_fig1 <- heat_tree_matrix(filt_mod_pseq_obj_with_parents,
                              data = "compare_tax_abund",
                              node_size = n_obs,
                              node_label = taxon_names,
                              node_color = log2_median_ratio,
                              node_color_range = diverging_palette(),
                              node_color_trans = "area",
                              node_color_interval = range(filt_mod_pseq_obj_with_parents$data$compare_tax_abund$log2_median_ratio, finite = TRUE),
                              node_size_axis_label = "Taxa counts",
                              node_color_axis_label = "Log2 median ratio",
                              layout = "davidson-harel",
                              key_size = 0.6,
                              seed = 314,
                              aspect_ratio = 1,
                              initial_layout = "reingold-tilford",
                              repel_labels=T,
                              output_file = "pseq_obj_filt_node_color_trans_linear.pdf",
                              verbose = TRUE)

filt_mod_pseq_with_parents_fig1

#save.image("heat-tree-Lee.RData")

