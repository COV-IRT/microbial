
  ## this is working off of the "dmm modeling troubleshooting" file posted here (https://covirt.slack.com/archives/D011KGN7293/p1607119487014500)
  # base is your code, but most comments will be mine in this doc

library(tidyverse)
library(phyloseq)

raw<-as_tibble(read.table("Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
dim(raw)

colnames(raw)<-gsub("NA_tax","unclass", colnames(raw)) %>% str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

### it seems something wasn't communicated clearly :/ sorry!
    ## as currently written, this seems to be keeping "term", "arc", "bac", "euk", "vir", and "unclassified" for each, where "term" itself holds all of them
    # here's the code as it was, and then breaking it down afterwards

df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts")) %>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","termteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)") %>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)

  # breaking down a subset to see it
raw_sub <- raw %>% select(1:28)
  # this holds 2 samples, CRR119894 and CRR119895
colnames(raw_sub)

  # this part drops the percent columns
raw_df_p1 <- raw_sub %>% select(GO_term, namespace, depth, name, ends_with("_counts"))
  # this part pivots it out
raw_df_p2 <- raw_df_p1 %>% pivot_longer(cols = -c(GO_term, namespace, depth, name), names_to = c("sample", "type", "abund"), names_pattern = "(.*)_(.*)_(.*)")
  # this part drops the "abund" column, and filters out based on having a count > 1
raw_df_p3 <- raw_df_p2 %>% select(-abund) %>% filter(value > 1)
  # this part converts it back to wie format, and adds zeroes for those needed
raw_df_p4 <- raw_df_p3 %>% pivot_wider(names_from = sample, values_from = value, values_fill = 0)

  # here we are back to our regular format, but the total for each sample should match what the total was for just the "term" column for that sample
raw_df_p4 %>% select(CRR119894, CRR119895) %>% colSums()
    # CRR119894 CRR119895
    # 125394264  49260762
raw_sub %>% select(CRR119894_term_counts, CRR119895_term_counts) %>% colSums()
    # CRR119894_term_counts CRR119895_term_counts
    #              62699648              24631944

  # almost exactly half of what happened above (e.g. the counts for "term", "bacteria", "euk", etc. were all being added before). I'm sure it's not exactly doubled/half because those with 1 count were removed
  # sorry!! it seems this wasn't communicated clearly, the "...term_counts" columns are the total GO term counts for each sample, the taxa ones that follow are all a portion of those total counts

  ## Starting over, reading in fresh
raw <- as_tibble(read.table("Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))

  # if we want a table of just all counts, here's one way
term_counts_df <- raw %>%
  select(GO_term, namespace, depth, name, ends_with("_term_counts"))

  # dropping the "_term_counts" part of the headers
colnames(term_counts_df) <- gsub("_term_counts", "", colnames(term_counts_df))
     # and swapping the _ for . in that one (i'm guessing to match a metadata file, BUT the "Combined_BALF_GO_Terms_metadata2.txt" metadata file read in below has it as an "_" too)
colnames(term_counts_df) <- gsub("NC1_SRR7796663", "NC1.SRR7796663", colnames(term_counts_df))

  # adding a column for rowsum of sample counts so can filter like done above
term_counts_df_with_row_sums <- term_counts_df
term_counts_df_with_row_sums$sum_across_samples <- term_counts_df %>% select(c(5:ncol(term_counts_df))) %>% rowSums()

  # can filter here if we want, like only keeping those > 1 count in any sample
term_counts_df_filtered <- term_counts_df_with_row_sums %>% filter(sum_across_samples > 1) %>% select(- sum_across_samples)

dim(term_counts_df) # 47,550
dim(term_counts_df_filtered) # 27,653

  # if filtering and keeping only those with at least 1 count, it's only a few hundred more
term_counts_df_filtered <- term_counts_df_with_row_sums %>% filter(sum_across_samples > 0) %>% select(- sum_across_samples)
dim(term_counts_df_filtered) # 28,169

  # copying to object 'df' so can use your code below
df <- term_counts_df_filtered

df$depth<-as.character(df$depth)
#SIDE NOTE:There are multiple processes and values for a single sample so you cant convert the sample to columns
  #### NOTE FROM MIKE LEE: not sure what this note is about, as the samples are columns currently, but maybe it is different now that we aren't retaining multiple for each sample (i.e. no longer keeping the total counts AND taxa breakdown counts, as before)

  ## splitting these after some changes made below (dropping a sample that isn't in the metadata currently)
# bio<-filter(df, namespace=="biological_process")
# mol<-filter(df, namespace=="molecular_function")


#REMEMBER TO FIX THIS AGAIN LATER 4 DEC 2020
# TROUBLESHOOTING CODE ONLY
#######################################################

####There are multiple processes and values for a single sample so you cant convert the sample to columns
#Old WAY

# #make individual tibbles for each type (bac, euk, term, arc, vir, etc)
# bio_term<-bio%>%filter(type=="bac")%>%select(-type)
# mol_term<-mol%>%filter(type=="bac")%>%select(-type)
#
# #subselect tibbles for only the counts and go terminology
# bio_term_counts<-bio_term%>%select(-c(namespace,depth,name))
# bio_term_tax<-bio_term%>%select(GO_term,namespace,depth,name)
# mol_term_counts<-mol_term%>%select(-c(namespace,depth,name))
# mol_term_tax<-mol_term%>%select(GO_term,namespace,depth,name)
#
# #convert them to dataframes for downstream import to phylsoeq
# bio_term_counts<-data.frame(bio_term_counts, row.names=1)
# bio_term_tax<-data.frame(bio_term_tax, row.names=1)
# mol_term_counts<-data.frame(mol_term_counts, row.names=1)
# mol_term_tax<-data.frame(mol_term_tax, row.names=1)
#
# #convert the dataframes into phyloseq formats
# bio_term_counts_phy <- otu_table(bio_term_counts, taxa_are_rows=TRUE)
# bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax), errorIfNULL=TRUE)
# mol_term_counts_phy<-otu_table(mol_term_counts, taxa_are_rows = T)
# mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax), errorIfNULL = T)

####################################################
# COMMENTED OUT ORIGINAL CODE
#############################################
# bio_term<-bio%>%filter(type%in%c("bac","arc","vir"))%>%
#   select(-type) %>%
#   group_by(GO_term,namespace,depth,name)%>%
#   summarise(across(.cols = where(is.numeric), sum))
#
# mol_term<-mol%>%
#   filter(type%in%c("bac","arc","vir"))%>%
#   select(-type)%>%
#   group_by(GO_term,namespace,depth,name)%>%
#   summarise(across(.cols = where(is.numeric), sum))
# bio_term_counts<-bio_term%>%
#   select(-c(namespace,depth,name))
# bio_term_tax<-bio_term%>%
#   select(GO_term,namespace,depth,name)
# mol_term_counts<-mol_term%>%
#   select(-c(namespace,depth,name))
# mol_term_tax<-mol_term%>%
#   select(GO_term,namespace,depth,name)

# bio_term_counts$namespace<-NULL
# bio_term_counts$depth<-NULL
# mol_term_counts$namespace<-NULL
# mol_term_counts$depth<-NULL

# bio_term_counts_df<-data.frame(bio_term_counts, row.names=1)
# bio_term_tax_df<-data.frame(bio_term_tax, row.names=1)
# mol_term_counts_df<-data.frame(mol_term_counts, row.names=1)
# mol_term_tax_df<-data.frame(mol_term_tax, row.names=1)
#
# dim(bio_term_counts)
# dim(bio_term_tax)
# dim(mol_term_counts)
# dim(mol_term_tax)
# #mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax), errorIfNULL = T)
#
# bio_term_counts_phy <- otu_table(bio_term_counts_df, taxa_are_rows=TRUE)
# bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax_df), errorIfNULL=TRUE)
# mol_term_counts_phy<-otu_table(mol_term_counts_df, taxa_are_rows = T)
# mol_term_tax_phy<-tax_table(as.matrix(mol_term_tax_df), errorIfNULL = T)
#################################################################


bio_term_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata2.txt", header = T, sep = "\t", row.names = 1))
dim(bio_term_sam) # 167 70
    # this says there are 167 samples, but we have 168 in our df above
  # seeing what's missing and/or different
samples_from_count_tab <- colnames(df)[5:ncol(df)]
length(samples_from_count_tab) # 168

samples_from_metadata_tab <- row.names(bio_term_sam)
length(samples_from_metadata_tab) # 167

setdiff(samples_from_count_tab, samples_from_metadata_tab)
  # "NC1.SRR7796663" "SRR10571760"
    # in the metadata table this NC1 has the underscore, so changing that back here


colnames(df) <- gsub("NC1.SRR7796663", "NC1_SRR7796663", colnames(df))
samples_from_count_tab <- colnames(df)[5:ncol(df)]
length(samples_from_count_tab) # 168

samples_from_metadata_tab <- row.names(bio_term_sam)
length(samples_from_metadata_tab) # 167
setdiff(samples_from_count_tab, samples_from_metadata_tab)
  # "SRR10571760" is missing from the metadata table

  ## for now, dropping it from counts table as I don't know what info is needed or not for later processing
df_SRR10571760_dropped <- df %>% select(- SRR10571760)
    # checking again all are the same now
samples_from_count_tab <- colnames(df_SRR10571760_dropped)[5:ncol(df_SRR10571760_dropped)]
length(samples_from_count_tab) # 167

samples_from_metadata_tab <- row.names(bio_term_sam)
length(samples_from_metadata_tab) # 167

setdiff(samples_from_count_tab, samples_from_metadata_tab)

 # checking a bit more robustly, this would come back false if there is a sample name that doesn't have a match in the other
all.equal(sort(samples_from_count_tab), sort(samples_from_metadata_tab))
   # TRUE

  # ah, i see you change this one here, now i'm guessing this was due to how you were parsing the sample names based on the "_" as a delimiter above (which i didn't end up doing, so skipping this)
# rownames(bio_term_sam)<-rownames(bio_term_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")


bio<-filter(df_SRR10571760_dropped, namespace=="biological_process")
mol<-filter(df_SRR10571760_dropped, namespace=="molecular_function")

bio_term_sam$accession<-rownames(bio_term_sam)

bio_term_sam$outcome<-bio_term_sam$outcome%>%
  str_replace_all("recovered", "Recovered")%>%
  str_replace_all("deceased","Deceased")%>%
  str_replace_all('stabilized',"Stabilized")

###############################################

  # this sex column when we read it in has "na" and "NA" in it, which when running this code below, is then leaving some as the special R object NA, and some as "<NA>" as a string:
# bio_term_sam$sex<-bio_term_sam$sex%>%
#   str_replace_all("M", "male")%>%
#   str_replace_all("F", "female")%>%
#   str_replace_all("na", "<NA>")
  # we can see both "NA" and "<NA>" values in there, and there are 3 factor levels (one including the string "na"), e.g.
bio_term_sam$sex
# [1] male   female male   male   female female female female male   female female female F      <NA>   M      M      M      F      <NA>   M      F      F      F      M
# [25] M      M      M      <NA>   <NA>   <NA>   <NA>   M      M      F      F      M      F      female male   male   male   female male   male   male   female female male
# [49] male   male   male   male   male   male   male   female female male   male   male   male   male   male   male   male   female female female female male   female male
# [73] male   male   male   female female male   female male   male   male   male   female female <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>
# [97] <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   female male   female male   female male   female male   male   female male   female <NA>   <NA>   <NA>
# [121] <NA>   na     female female female male   female female male   female female female male   female female male   female female male   female female female male   male
# [145] male   female na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na     na
# Levels: F female M male na

  # so better to read it in and tell it there are other NA characters we want to be treated as special NAs
bio_term_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata2.txt", header = T, sep = "\t", row.names = 1, na.strings=c("NA", "na")))

  # now if we look at that column, there are no mixed ones in there:
bio_term_sam$sex
# [1] male   female male   male   female female female female male   female female female F      <NA>   M      M      M      F      <NA>   M      F      F      F      M
# [25] M      M      M      <NA>   <NA>   <NA>   <NA>   M      M      F      F      M      F      female male   male   male   female male   male   male   female female male
# [49] male   male   male   male   male   male   male   female female male   male   male   male   male   male   male   male   female female female female male   female male
# [73] male   male   male   female female male   female male   male   male   male   female female <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>
# [97] <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   female male   female male   female male   female male   male   female male   female <NA>   <NA>   <NA>
# [121] <NA>   <NA>   female female female male   female female male   female female female male   female female male   female female male   female female female male   male
# [145] male   female <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>   <NA>
# Levels: F female M male

  # here still changing all male/female to be the same
bio_term_sam$sex <- bio_term_sam$sex %>%
  str_replace_all("M", "male") %>%
  str_replace_all("F", "female")

  # and still adding the accession column you had added
bio_term_sam$accession<-rownames(bio_term_sam)


### making phyloseq objects sort of based on how you were doing them in above commented-out code, but also some other steps i missed or didn't see probably because i'm jumping around, and then taking off again with your code
bio_term_counts <- bio %>% select(-c(namespace, depth, name)) # in what i saw above, these didn't keep the GO_term IDs, which i think need to be the rownames here when we make the phyloseq object, so they match up with the tax object you're making (so it's able to link the counts to the "tax" info)
bio_term_counts_df <- bio_term_counts %>% column_to_rownames("GO_term")

bio_term_tax <- bio %>% select(GO_term, namespace, depth, name)
bio_term_tax_df <- bio_term_tax %>% column_to_rownames("GO_term")


mol_term_counts <- mol %>% select(-c(namespace, depth, name))
mol_term_counts_df <- mol_term_counts %>% column_to_rownames("GO_term")

mol_term_tax <- mol %>% select(GO_term, namespace, depth, name)
mol_term_tax_df <- mol_term_tax %>% column_to_rownames("GO_term")



bio_term_counts_phy <- otu_table(bio_term_counts_df, taxa_are_rows=TRUE)
bio_term_tax_phy <- tax_table(as.matrix(bio_term_tax_df), errorIfNULL=TRUE)

mol_term_counts_phy <- otu_table(mol_term_counts_df, taxa_are_rows=TRUE)
mol_term_tax_phy <- tax_table(as.matrix(mol_term_tax_df), errorIfNULL=TRUE)

bio_term_pseq <- phyloseq(bio_term_counts_phy, bio_term_tax_phy, sample_data(bio_term_sam))
mol_term_pseq<-phyloseq(mol_term_counts_phy,mol_term_tax_phy, sample_data(bio_term_sam))


#### i'm not sure why we'd be merging these phyloseq objects here, and if we want them together, why we didn't just *not* split them above
#### but running through as similarly as possible just to hopefully find the same problem you hit

term_pseq<-merge_phyloseq(bio_term_pseq,mol_term_pseq)
term_pseq# [ 14581 taxa and 167 samples ]  ## # got 28169 taxa and 167 samples when Mike Lee got here from doing the above (i think would could be related to them not going in with GO ID as row names)

  # i'm skipping this for now as not sure why we'd be filtering it out
# filtme<-c("GO:0003674")
# term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
# filtme<-c("GO:0008150")
# term_pseq <- prune_taxa(taxa=taxa_names(term_pseq)!=filtme, term_pseq)
# term_pseq

term_pseq_no_neg<-subset_samples(term_pseq, sample_type!="neg_control")
term_pseq_no_neg# [ 14579 taxa and 162 samples ] ## # 28169 and 162 when Mike Lee got here
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, sample_type!="Unknown")
term_pseq_no_neg#  [ 14579 taxa and 141 samples ]: ## # 28169 and 141 for ML
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, case!="Control_Sick")
term_pseq_no_neg# [ 14597 taxa and 105 samples ] ## # 28169 and 105 for ML
term_pseq_no_neg<-subset_samples(term_pseq_no_neg,publication!="Michalovich")
term_pseq_no_neg# [ 14597 taxa and 102 samples ] ## # 28169 and 102 for ML
term_pseq_no_neg<-subset_samples(term_pseq_no_neg, bioproject!="PRJNA605907")
term_pseq_no_neg# [ 14597 taxa and 86 samples ] ## # 28169 and 86 for ML
term_pseq_no_neg<-prune_taxa(taxa = taxa_sums(term_pseq_no_neg)>0, x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ] ## # 27526 and 86 for ML
term_pseq_no_neg<-prune_samples(samples = sample_sums(term_pseq_no_neg)>0,x = term_pseq_no_neg)
term_pseq_no_neg# [ 13534 taxa and 86 samples ] ## # 27526 and 86 for ML


# BiocManager::install("microbiome") # ML install
packageVersion("microbiome") # 1.6.0 ML

pseq <- term_pseq_no_neg

##### i'm not sure i'm looking at the right stuff below, so just trying to plug this into DMM to see if it errors
# some of what i changed above might affect the problem you've been hitting
dat <- abundances(pseq)
count <- as.matrix(t(dat))

# BiocManager::install("DirichletMultinomial") # ML install
packageVersion("DirichletMultinomial") # 1.26.0 ML
library(DirichletMultinomial)

  # just seeing if this errors out, not intending to let run
fit <- lapply(1:3, dmn, count = count, verbose = TRUE)
    # this starts and prints out "Soft kmeans" and then "Expectation Maximization setup"
       # would there be an error immediately?


##### ML didn't go any further than here #####

#
# # To speed up, only consider the core taxa
# # that are prevalent at 0.1% relative abundance in 50% of the samples
# # (note that this is not strictly correct as information is
# # being discarded; one alternative would be to aggregate rare taxa)
# pseq.comp <- microbiome::transform(pseq, "compositional")
# taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
# pseq <- prune_taxa(taxa, pseq)
# pseq #[ 133 taxa and 86 samples ]
# # Pick the OTU count matrix
# # and convert it into samples x taxa format
# dat <- abundances(pseq)
# count <- as.matrix(t(dat))
# rowSums(count)
# colSums(count)
# #Fit the DMM model. Let us set the maximum allowed number of community types to 3 to speed up the example.
# #
# #install.packages('tictoc')
#
#
# library(DirichletMultinomial)
# library(tictoc)
# tic("dmm_modeling")
# print("dmm modeling starting")
# fit <- lapply(1:3, dmn, count = count, verbose=TRUE)
# print("dmm modeling finished")
# toc()
#
# #The is strangely taking a very long time for a physeq object that only has
# #133 terms and 86 samples.  Perhaps the model is not limited by the number of terms, as much as it is
# #limited by the size and disparity of the counts in the dataset   its 11:45
#
# #Check model fit with different number of mixture components using standard information criteria
# lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
# aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
# bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
# plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
# #Pick the optimal model
# best <- fit[[which.min(unlist(lplc))]]
# #Mixture parameters pi and theta
# mixturewt(best)
# #Sample-component assignments
# ass <- apply(mixture(best), 1, which.max)
# # Contribution of each taxonomic group to each component
# for (k in seq(ncol(fitted(best)))) {
#   d <- melt(fitted(best))
#   colnames(d) <- c("OTU", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange OTUs by assignment strength
#     arrange(value) %>%
#     mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
#     # Only show the most important drivers
#     filter(abs(value) > quantile(abs(value), 0.8))
#
#   p <- ggplot(d, aes(x = OTU, y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     labs(title = paste("Top drivers: community type", k))
#   print(p)
# }