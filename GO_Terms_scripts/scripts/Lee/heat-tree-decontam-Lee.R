library(tidyverse)
library(phyloseq)
# setwd("c:/github/microbial/step2_kraken2_analysis/")

# if wanting to load the objects
load("heat-tree-Jochum.RData")

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


pseq<-phyloseq(otu_table(bac_tax_counts_tab,taxa_are_rows = T),tax_table(as.matrix(bac_tax_info_tab)),sample_data(target_sample_info_tab))
pseq

### diverged from your initial code in some places starting here at the decontam part
library(decontam)
#export the shen et al. samples for parsing with decontam
pseq_shen<-subset_samples(pseq,bioproject=="CRA002476")
pseq_not_shen<-subset_samples(pseq,bioproject!="CRA002476")
pseq_shen

# getting count table of just shen
count_df_shen<-abundances(pseq_shen)

# getting sample info for just shen
sam_shen<-meta(pseq_shen)

# getting T/F vector of which are negative controls
vector_for_decontam <- sam_shen$case=="Control_Neg"
vector_for_decontam

# running decontam
contam_df_shen <- isContaminant(t(count_df_shen), neg=vector_for_decontam)
table(contam_df_shen$contaminant) # this gives me 4,527 as contaminants, previous value in following parentheses (identified 4644 as contaminants)

# getting vector holding the identified contaminant IDs
contam_shen <- row.names(contam_df_shen[contam_df_shen$contaminant == TRUE, ])

count_df_shen<-as_tibble(count_df_shen,rownames = "taxid")
count_df_shen

clean_count_df_shen<-count_df_shen%>%filter(!taxid%in%contam_shen)
dim(clean_count_df_shen)

a<-clean_count_df_shen%>%select(-taxid) # though my value of contaminants above was different, this object "a" had the same count before i altered it (4,527), so maybe that was just an old note above
b<-count_df_shen%>%select(-taxid)
perc_remain<-(colSums(a)/colSums(b))*100

write.table(perc_remain,"percent_remaining_post_decontam.tsv",sep = "\t",quote = F,row.names = F)
tax<-as.tibble(bac_tax_info_tab,rownames="taxid")
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

###### Whoops, i think something funny happened here, it seems we're taking the contaminants and sticking them into pseq_shen_decontam, which is pretty confusing, haha
pseq_shen_decontam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)
pseq_shen_decontam
###### Then, we are merging the shen contaminants with the rest of the data
pseq_decontam<-merge_phyloseq(pseq_not_shen,pseq_shen_decontam)
###### This object was used moving forward and so i'm pretty sure the final tree before was made a little wonky with regard to the shen samples
###### here we'll switch this, so we are keeping the non-contaminants in the pseq_shen_decontam object
pseq_shen_decontam<-prune_taxa(taxa = not_contaminants$taxid, pseq_shen)
pseq_shen_decontam # now it holds 4,134
pseq_decontam<-merge_phyloseq(pseq_not_shen, pseq_shen_decontam)
    ### moving forward with this updated one, should have 4,134


##########
  # out of curiousity, checking how much is retained if we remove those Shen-contaminant-flagged taxids from all samples
shen_decontamd_counts_tab <- otu_table(pseq_decontam) %>% data.frame()
shen_contaminants_removed_from_all_counts_tab <- shen_decontamd_counts_tab[!row.names(shen_decontamd_counts_tab) %in% contam_shen, ]
dim(shen_decontamd_counts_tab)
dim(shen_contaminants_removed_from_all_counts_tab)
colSums(shen_contaminants_removed_from_all_counts_tab) / colSums(shen_decontamd_counts_tab) * 100
##########

pseq_decontam_no_neg<-subset_samples(physeq = pseq_decontam,case!="Control_Neg")

library(microbiome)
pseq_decontam_no_neg_core<-core(x = pseq_decontam_no_neg,detection = 0,prevalence = 50/100)
pseq_decontam_no_neg_core # previously [ 566 taxa and 86 samples ]; now 683 taxa, 86 samples




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


head(tab_for_taxa_object)
head(hmp_otus)[, 1:10]
str(tab_for_taxa_object)
   # i think all's ok, let's give it a shot

pseq_obj<-parse_phyloseq(pseq_decontam_no_neg_core)
pseq_obj
taxa_obj <- parse_tax_data(tab_for_taxa_object,
                           class_cols = "lineage",
                           class_sep = ";",
                           class_regex = "^([c-s]{1})__(.*)$",
                           class_key = c(tax_rank = "info", tax_name = "taxon_name"))

taxa_obj
# looks ok so far, following along with same example page to filter low abundance things, just picking anything

taxa_obj_filt <- taxa::filter_taxa(taxa_obj, n_obs > 5)
pseq_obj_filt <- taxa::filter_taxa(pseq_obj, n_obs > 5)

taxa_obj_filt
pseq_obj_filt
# this is where we normalize across samples (i'm curious what this spits out when there are negative values in the table like from the vst, ha)
#################################################################################
pseq_obj_filt$data$tax_proportions <- calc_obs_props(pseq_obj_filt, "otu_table") # call this by the otu-table

taxa_obj_filt$data$tax_proportions <- calc_obs_props(taxa_obj_filt, "tax_data")
########################################################################################

# i don't yet understand why the calc_taxon_abund() step is dropping us down to 755, when all of our taxa are unique already... but oh well

pseq_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_obj_filt, "tax_proportions", cols = target_samples_no_neg)
taxa_obj_filt$data$tax_abund <- calc_taxon_abund(taxa_obj_filt, "tax_proportions", cols = target_samples)

# not sure if we want the groups argument here or not
pseq_obj_filt$data$tax_occ <- calc_n_samples(pseq_obj_filt, "tax_abund")
taxa_obj_filt$data$tax_occ <- calc_n_samples(taxa_obj_filt, "tax_abund")

# doing the compare now with the abund table
pseq_obj_filt$data$compare_tax_abund <- compare_groups(pseq_obj_filt,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_obj_filt$data$sample_data$sample_type)



taxa_obj_filt$data$compare_tax_abund <- compare_groups(taxa_obj_filt,
                                                             "tax_abund",
                                                             cols = target_samples,
                                                             groups = target_sample_info_tab$sample_type)


range(pseq_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# looking at -5.2 to 4.2

range(taxa_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# we're looking at -6.8 to 7.3

# trying a plot (commenting out node_label = taxon_names makes them plot much faster once saved to an object, good for when trying different things)

fig1 <- heat_tree_matrix(taxa_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "linear",
                         node_color_interval = c(-7.5,7.5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "node_color_trans_linear.pdf",
                         verbose = TRUE)
pseq_fig1 <- heat_tree_matrix(pseq_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "linear",
                         node_color_interval = c(-4,4),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "node_color_trans_linear.pdf",
                         verbose = TRUE)

pseq_fig1
fig2 <- heat_tree_matrix(taxa_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-7.5,7.5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "node_color_trans_area.pdf",
                         verbose = TRUE)

pseq_fig2 <- heat_tree_matrix(pseq_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-4,4),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "node_color_trans_area.pdf",
                         verbose = TRUE)

pseq_fig2
## it looks to me like the default "area" node_color_trans gives a better pop than the "linear"

## i think we might want to slim it down more, and i think we can do that much more easily outside of metacoder, just like
## we made the object above with what we wanted, we can trim them down first and then not have to worry about trimming them down in there (which is confusing me, ha)
## although maybe we can do it conveniently enough with the supertaxa() function: https://github.com/ropensci/taxa#supertaxa
## i haven't tried it yet, but it seems that may help. Looks like we can pick a rank and tell it to bump everything up to that level

## trying the supertaxa() function prior to prevalence filtering
## first just making the fuller ones
pseq_decontam_no_neg

## creating taxa object
pseq_full_taxa_obj <- parse_phyloseq(pseq_decontam_no_neg) # i wonder if we should set the NAs to blank strings prior to this, given the warning message about NAs being in there and how their example has nothing instead
pseq_full_taxa_obj

pseq_full_taxa_obj_filt <- taxa::filter_taxa(pseq_full_taxa_obj, n_obs > 5) # this took us down to 731 taxa
pseq_full_taxa_obj_filt

pseq_full_taxa_obj_filt$data$tax_proportions <- calc_obs_props(pseq_full_taxa_obj_filt, "otu_table")

pseq_full_taxa_obj_filt$data$tax_abund <- calc_taxon_abund(pseq_full_taxa_obj_filt, "tax_proportions", cols = target_samples_no_neg)

# not sure if we want the groups argument here or not
pseq_full_taxa_obj_filt$data$tax_occ <- calc_n_samples(pseq_full_taxa_obj_filt, "tax_abund")

# doing the compare now with the abund table
pseq_full_taxa_obj_filt$data$compare_tax_abund <- compare_groups(pseq_full_taxa_obj_filt,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_full_taxa_obj_filt$data$sample_data$sample_type)

range(pseq_full_taxa_obj_filt$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# looking at -6.4 to 9.6

pseq_full_taxa_obj_fig1 <- heat_tree_matrix(pseq_full_taxa_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
#                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "linear",
                         node_color_interval = c(-6.5,9.5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_node_color_trans_linear.pdf",
                         verbose = TRUE)


pseq_full_taxa_obj_fig2 <- heat_tree_matrix(pseq_full_taxa_obj_filt,
                         data = "compare_tax_abund",
                         node_size = n_obs,
#                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-6.5,9.5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_node_color_trans_area.pdf",
                         verbose = TRUE)

# those are pretty full of course

####### trying supertaxa business finally, i imagine we should do it before adding all the things
pseq_full_taxa_obj_for_supertaxa_try <- parse_phyloseq(pseq_decontam_no_neg) # i wonder if we should set the NAs to blank strings prior to this, given the warning message about NAs being in there and how their example has nothing instead
pseq_full_taxa_obj_for_supertaxa_try

pseq_full_taxa_obj_supertaxa <- supertaxa(pseq_full_taxa_obj_for_supertaxa_try)
     ### okay, that's not doing what I though it might, but i realize i think we can do this with the taxa::filter_taxa function too
pseq_full_taxa_obj_genus_collapsed <- taxa::filter_taxa(pseq_full_taxa_obj_for_supertaxa_try, taxon_ranks == "genus", supertaxa = TRUE)
pseq_full_taxa_obj_genus_collapsed # genus by itself is still 2,010, so let's try higher
pseq_full_taxa_obj_family_collapsed <- taxa::filter_taxa(pseq_full_taxa_obj_for_supertaxa_try, taxon_ranks == "family", supertaxa = TRUE)
pseq_full_taxa_obj_family_collapsed # this droped us to 646, going one more
pseq_full_taxa_obj_order_collapsed <- taxa::filter_taxa(pseq_full_taxa_obj_for_supertaxa_try, taxon_ranks == "order", supertaxa = TRUE)
pseq_full_taxa_obj_order_collapsed # this dropped us to 271, let's see how possible that is

pseq_full_taxa_obj_order_collapsed$data$tax_proportions <- calc_obs_props(pseq_full_taxa_obj_order_collapsed, "otu_table")

pseq_full_taxa_obj_order_collapsed$data$tax_abund <- calc_taxon_abund(pseq_full_taxa_obj_order_collapsed, "tax_proportions", cols = target_samples_no_neg)

# not sure if we want the groups argument here or not
pseq_full_taxa_obj_order_collapsed$data$tax_occ <- calc_n_samples(pseq_full_taxa_obj_order_collapsed, "tax_abund")

# doing the compare now with the abund table
pseq_full_taxa_obj_order_collapsed$data$compare_tax_abund <- compare_groups(pseq_full_taxa_obj_order_collapsed,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_full_taxa_obj_order_collapsed$data$sample_data$sample_type)

range(pseq_full_taxa_obj_order_collapsed$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# looking at -6.4 to 9.6


pseq_full_taxa_obj_order_collapsed_fig2 <- heat_tree_matrix(pseq_full_taxa_obj_order_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-6.5,9.5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_order_node_color_trans_area.pdf",
                         verbose = TRUE)

### not bad, let's look at one higher just for the f of it
pseq_full_taxa_obj_class_collapsed <- taxa::filter_taxa(pseq_full_taxa_obj_for_supertaxa_try, taxon_ranks == "class", supertaxa = TRUE)
pseq_full_taxa_obj_class_collapsed # this dropped us to 107, let's see how possible that is

pseq_full_taxa_obj_class_collapsed$data$tax_proportions <- calc_obs_props(pseq_full_taxa_obj_class_collapsed, "otu_table")

pseq_full_taxa_obj_class_collapsed$data$tax_abund <- calc_taxon_abund(pseq_full_taxa_obj_class_collapsed, "tax_proportions", cols = target_samples_no_neg)

# not sure if we want the groups argument here or not
pseq_full_taxa_obj_class_collapsed$data$tax_occ <- calc_n_samples(pseq_full_taxa_obj_class_collapsed, "tax_abund")

# doing the compare now with the abund table
pseq_full_taxa_obj_class_collapsed$data$compare_tax_abund <- compare_groups(pseq_full_taxa_obj_class_collapsed,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_full_taxa_obj_class_collapsed$data$sample_data$sample_type)

range(pseq_full_taxa_obj_class_collapsed$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# looking at -4 to 7


pseq_full_taxa_obj_class_collapsed_fig <- heat_tree_matrix(pseq_full_taxa_obj_class_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-4,7),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_class_node_color_trans_area.pdf",
                         verbose = TRUE)

  # i don't think n_obs for node size is that informative for us here (because our things are already unique-d, as opposed to the case when giving this taxonomy of amplicon data),
  # going to try number of samples (in the tax_occ table in there)
pseq_full_taxa_obj_class_collapsed_fig <- heat_tree_matrix(pseq_full_taxa_obj_class_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_samples,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-4,7),
                         node_size_axis_label = "Number of samples detected",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_class_node_color_trans_area_node_size_n_samples.pdf",
                         verbose = TRUE)

  # hmm, that could be scaled better, let's try
pseq_full_taxa_obj_class_collapsed_fig <- heat_tree_matrix(pseq_full_taxa_obj_class_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_samples,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_size_trans = "log10",
                         node_color_interval = c(-4,7),
                         node_size_axis_label = "Number of samples detected",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "full_class_node_color_trans_area_node_size_n_samples_scaled.pdf",
                         verbose = TRUE)

   ## guess not, i think it's just cause we bumped up so many ranks, maybe filtering is the way to go first like you were doing, let's try that

sample_data(pseq_decontam_no_neg)$case %>% table
  # our smallest group has 25, so let's see what happens requiring things to be in at least that many

new_pseq_decontam_no_neg_core<-core(x = pseq_decontam_no_neg, detection = 0, prevalence = 25/86)
  # that cut down to 1,858

new_pseq_decontam_no_neg_core_taxa_obj <- parse_phyloseq(new_pseq_decontam_no_neg_core)
new_pseq_decontam_no_neg_core_taxa_obj

pseq_taxa_obj_genus_collapsed <- taxa::filter_taxa(new_pseq_decontam_no_neg_core_taxa_obj, taxon_ranks == "genus", supertaxa = TRUE)
pseq_taxa_obj_genus_collapsed # genus by itself is still 1,006, so let's try higher
pseq_taxa_obj_family_collapsed <- taxa::filter_taxa(new_pseq_decontam_no_neg_core_taxa_obj, taxon_ranks == "family", supertaxa = TRUE)
pseq_taxa_obj_family_collapsed # 405, now we're gettin somewhere, let's do order
pseq_taxa_obj_order_collapsed <- taxa::filter_taxa(new_pseq_decontam_no_neg_core_taxa_obj, taxon_ranks == "order", supertaxa = TRUE)
pseq_taxa_obj_order_collapsed # 172, with some clean-up, this might work well

pseq_taxa_obj_order_collapsed$data$tax_proportions <- calc_obs_props(pseq_taxa_obj_order_collapsed, "otu_table")

pseq_taxa_obj_order_collapsed$data$tax_abund <- calc_taxon_abund(pseq_taxa_obj_order_collapsed, "tax_proportions", cols = target_samples_no_neg)

pseq_taxa_obj_order_collapsed$data$tax_occ <- calc_n_samples(pseq_taxa_obj_order_collapsed, "tax_abund")

# doing the compare now with the abund table
pseq_taxa_obj_order_collapsed$data$compare_tax_abund <- compare_groups(pseq_taxa_obj_order_collapsed,
                                                       "tax_abund",
                                                       cols = target_samples_no_neg,
                                                       groups = pseq_taxa_obj_order_collapsed$data$sample_data$sample_type)

range(pseq_taxa_obj_order_collapsed$data$compare_tax_abund$log2_median_ratio, finite = TRUE)
# looking at -5 to 5


pseq_taxa_obj_order_collapsed_fig <- heat_tree_matrix(pseq_taxa_obj_order_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = taxon_names,
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-5,5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "subset_order_heattree.pdf",
                         verbose = TRUE)

  ### I don't hate that, let's try getting a vector of somes names we *don't* want to plot
  # first testing it works the way i think but saying not to plot "Bacteria"
test_fig <- heat_tree_matrix(pseq_taxa_obj_order_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = ifelse(taxon_names %in% "Bacteria", "", taxon_names),
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-5,5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "TEST-subset_order_heattree.pdf",
                         verbose = TRUE)

  # beautiful, let's see what we want to *not* plot
  # i don't want to remove anything below the phylum level for sure, and maybe the class level, but we'll see
  # let's look at the distribution of log2 diffs
compare_tab <- pseq_taxa_obj_order_collapsed$data$compare_tax_abund
  # removing infs
compare_tab <- compare_tab[is.finite(compare_tab$log2_median_ratio), ]
dim(compare_tab)
summary(compare_tab$log2_median_ratio)
hist(compare_tab$log2_median_ratio)
quantile(compare_tab$log2_median_ratio, seq(0,1,0.01))

sd(compare_tab$log2_median_ratio)

sum(compare_tab$log2_median_ratio > sd(compare_tab$log2_median_ratio)) # there are 66 above 1SD
sum(compare_tab$log2_median_ratio < -sd(compare_tab$log2_median_ratio)) # 36 below 1SD
   # and those will actually be fewer because some of them will be in there twice (due to the different contrasts in the table)
   # so let's try only displaying those (and the higher taxa ranks too)
   # due to wanting to keep the higher taxa ranks too, I think it'll be easier to get a list of what we want, and then invert it
taxon_ranks(pseq_taxa_obj_order_collapsed)
   # since this and the compare_tax_abund table are using the taxon IDs (aab, aad, etc.), going to work with them instead rather than converting all over the place (assuming we can pass those to the plot function - we'll figure that out then...)
   # so first getting all those down to the class rank
all_taxon_IDs_and_ranks <- taxon_ranks(pseq_taxa_obj_order_collapsed)
wanted_taxon_IDs_based_on_rank <- all_taxon_IDs_and_ranks[all_taxon_IDs_and_ranks %in% c("domain", "phylum", "class")] %>% names()
length(wanted_taxon_IDs_based_on_rank) # we're at 67 names already, that seems like a lot, may come back and try phylum level
   # now getting those above 1SD or below 1SD
wanted_taxon_IDs_based_on_diff_abund <- compare_tab[compare_tab$log2_median_ratio > sd(compare_tab$log2_median_ratio) | compare_tab$log2_median_ratio < -sd(compare_tab$log2_median_ratio), "taxon_id"] %>% pull() %>% unique()
length(wanted_taxon_IDs_based_on_diff_abund) # 65, lets see how many together

all_wanted_taxon_IDs <- union(wanted_taxon_IDs_based_on_rank, wanted_taxon_IDs_based_on_diff_abund)
length(all_wanted_taxon_IDs) # 112, let's see how crowded it looks, but i think we'll be slimming down some more
  # inverting now to get the vector we'll pass to the plotting call
unwanted_taxon_IDs <- names(all_taxon_IDs_and_ranks)[!names(all_taxon_IDs_and_ranks) %in% all_wanted_taxon_IDs]

test_fig <- heat_tree_matrix(pseq_taxa_obj_order_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = ifelse(taxon_ids %in% unwanted_taxon_IDs, "", taxon_names),
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         node_color_trans = "area",
                         node_color_interval = c(-5,5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "TEST-subset_order_heattree.pdf",
                         verbose = TRUE)

## still too many too small, no respectable journal is going to let us get away with like 4-pt font on a figure..., going to up things to 2 SD
wanted_taxon_IDs_based_on_diff_abund <- compare_tab[compare_tab$log2_median_ratio > 2 * sd(compare_tab$log2_median_ratio) | compare_tab$log2_median_ratio < 2 * -sd(compare_tab$log2_median_ratio), "taxon_id"] %>% pull() %>% unique()
length(wanted_taxon_IDs_based_on_diff_abund) # 20, lets see how many together

all_wanted_taxon_IDs <- union(wanted_taxon_IDs_based_on_rank, wanted_taxon_IDs_based_on_diff_abund)
length(all_wanted_taxon_IDs) # 83
  # inverting now to get the vector we'll pass to the plotting call
unwanted_taxon_IDs <- names(all_taxon_IDs_and_ranks)[!names(all_taxon_IDs_and_ranks) %in% all_wanted_taxon_IDs]

pseq_taxa_obj_order_collapsed2<-filter_taxa(obj = pseq_taxa_obj_order_collapsed,wanted_taxon_IDs_based_on_diff_abund)
write.table(pseq_taxa_obj_order_collapsed2$data$tax_data,"diff_tax_data.tsv",sep = "\t")

test_fig <- heat_tree_matrix(pseq_taxa_obj_order_collapsed,
                         data = "compare_tax_abund",
                         node_size = n_obs,
                         node_label = ifelse(taxon_ids %in% unwanted_taxon_IDs, "", taxon_names),
                         node_color = log2_median_ratio,
                         node_color_range = diverging_palette(),
                         key_size = 0.7,
                         node_color_trans = "area",
                         node_color_interval = c(-5,5),
                         node_size_axis_label = "Taxa counts",
                         node_color_axis_label = "Log2 median ratio",
                         layout = "davidson-harel",
                         initial_layout = "reingold-tilford",
                         output_file = "TEST-subset_order_heattree.pdf",
                         verbose = TRUE)

test_fig
wanted_taxon_IDs_based_on_rank
## that's not terrible, I'm not going to get to the manual work in Affinity Designer on this figure at the moment,
   ## but i think we might want to make the map tree like a full page size, and then the 3 others being side-by-side beneath it or something
   ## they can be smaller and still do their job, but we really need the main map tree bigger, i will try that when i get to messing with fine-touches on the figure

   ## other thought to remember for now, I'm still not sold that n_obs for node size is telling us all that much, maybe we can chat about that next meeting to see if i'm just thinking of it wrong
   ## other, other thought, I'd like us to make sure at some point it's not doing anything weird with those NAs in the tax when coming from phyloseq,
      ## pulling out the data, and making it from scratch like i did following their example stuff just to make sure all the numbers and figures look the same in the end is one way


## last note for now while i'm remembering it, there are NAs still in the taxonomy, because there are stupid ones that are right in the middle, and not at the end where my little code above took care of them
## we can see them with this:
lineage_vector[grep("NA", x=lineage_vector)]
    ## some show up in the plot, and some nodes have no labels, not sure if that's another problem or not currently
    ## some of these are consistent enough it'd be easy to fix, like the Cyanobacteria all have NA for class, we can fix that and it cuts
    ## this list from 530 down to like 330, if there are others that are big chunks, maybe we could fix them all if we wanted



