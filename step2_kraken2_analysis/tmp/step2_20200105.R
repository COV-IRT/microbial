library(phyloseq)
library(tidyverse)
#setwd("/media/jochum00/Aagaard_Raid3/microbial/step2_kraken2_analysis/")
setwd("c:/github/microbial/step2_kraken2_analysis/")

#filter out only the 86 samples used in the GO term analysis
sam<-as_tibble(read.table("Combined_BALF_GO_Terms_metadata.tsv",header = T,sep = "\t"))
go_term_samples<-as_tibble(read.table("GO_TERM_samples.tsv",header=T,sep="\t"))
negs<-as_tibble(c("CRR125995",'CRR125996',"CRR125997","CRR125998"))
negs$accession<-negs$value
negs$value<-NULL
go_term_samples
keep<-bind_rows(go_term_samples,negs)
keep
sam<-left_join(keep,sam)
sam
#make a vector of the sample names
vec<-sam$accession
#filter out only the 86 samples used in the go term analysis
counts<-as_tibble(read.table("counts.tsv",header=T,sep="\t"))%>%
  pivot_longer(cols = -c(taxid),names_to="accession")%>%
  filter(accession %in% vec)%>%
  pivot_wider(names_from=accession)
#import the taxonomy and merge it up with the counts
tax<-as_tibble(read.table("taxonomy.tsv",header=T,sep="\t"))#,rownames = "taxid")
# convert to df
counts_df<-data.frame(counts,row.names = 1)
tax_df<-data.frame(tax,row.names = 1)
sam_df<-data.frame(sam,row.names=1)
#make a phyloseq obj
count_pseq<-otu_table(count_df,taxa_are_rows = T)
tax_pseq<-tax_table(as.matrix(tax_df),errorIfNULL = T)
sam_pseq<-sample_data(sam_df)
pseq<-phyloseq(count_pseq,tax_pseq,sam_pseq)
pseq<-subset_taxa(physeq = pseq,domain=="Bacteria")

########################################################################################
library(decontam)
#export the shen et al. samples for parsing with decontam
pseq_shen<-subset_samples(pseq,bioproject=="CRA002476")
pseq_not_shen<-subset_samples(pseq,bioproject!="CRA002476")

count_df_shen<-abundances(pseq_shen)
tax_table(pseq_shen) %>% as("matrix")# %>% subset(select = colSums(!is.na(.)) > 0) %>% as_tibble(rownames = "OTU") -> taxon.tbl

tax_df_shen<-as.matrix(tax_table(pseq_shen))
tax_df
tax_df_shen<-data.frame(tax_table(pseq_shen))
tax_df_shen
sam_shen<-meta(pseq_shen)
vector_for_decontam <- sam_shen$case=="Control_Neg"
vector_for_decontam                    
colSums(count_df_shen)
contam_df_shen <- isContaminant(t(count_df_shen), neg=vector_for_decontam)
table(contam_df_shen$contaminant) # identified 4644 as contaminants 

# getting vector holding the identified contaminant IDs
contam_shen <- row.names(contam_df_shen[contam_df_shen$contaminant == TRUE, ])
contam_shen<-as.factor(contam_shen)
contam_shen
tax$taxid<-as.factor(tax$taxid)
contaminants<-tax%>%filter(taxid %in% contam_shen)

write.table(x = contaminants,file = "Shen_contaminants.tsv",sep = "\t",row.names = F,quote = F)

contam_shen<-as.factor(contam_shen)
contaminants$taxid<-as.character(contaminants$taxid)
pseq_shen_decontam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)
pseq_shen_contam<-prune_taxa(taxa = contaminants$taxid,pseq_shen)

pseq_shen_decontam

pseq_decontam<-merge_phyloseq(pseq_not_shen,pseq_shen_decontam)
pseq_decontam_no_neg<-subset_samples(physeq = pseq_decontam,case!="Control_Neg")

pseq_decontam_no_neg

###############################################################################################
#vst normalization
library(DESeq2)

deseq_counts<-phyloseq_to_deseq2(pseq_decontam_no_neg, design = ~ NULL)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 






# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(pseq_decontam_no_neg)
vst_physeq <- phyloseq(vst_count_phy,tax_pseq, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="case") + 
  geom_point(size=1) + labs(col="case") + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") 





library(metacoder)

vst_physeq_bac<-subset_taxa(vst_physeq_bac,domain=="Bacteria")
vst_physeq_bac_prune<-core(x = vst_physeq_bac,detection = 0,prevalence = 50/100)
vst_physeq_bac_prune
min(taxa_sums(vst_physeq_bac_prune))
max(taxa_sums(vst_physeq_bac_prune))
hist(taxa_sums(vst_physeq_bac_prune))


vst_physeq_bac_prune
sam_no_neg<-sam%>%filter(case!="Control_Neg")
sam_no_neg
######################
obj<-parse_phyloseq(vst_physeq_bac_prune)

obj$data$otu_sqrt<-obj$data$otu_table%>%
  pivot_longer(-c(taxon_id,otu_id),names_to="accession",values_to="abundance")%>%
  mutate(sqrtabund=sqrt(abundance))%>%select(-abundance)%>%pivot_wider(names_from = accession,values_from=sqrtabund)

obj_filt<-filter_taxa(obj,n_obs>5) # remove taxa w/less than 5 observations for each taxon in a data set in a taxmap()
#object. This includes observations for the specific taxon and the observations of its subtaxa.obj_filt # 84 taxa

#For a given table in a taxmap object, convert one or more columns containing counts to proportions.
#This is meant to be used with counts associated with observations (e.g. OTUs), 
#as opposed to counts that have already been summed per taxon.
obj_filt$data$tax_abund <- calc_taxon_abund(obj_filt, "otu_table",
                                            groups = NULL,#sam_no_neg$sample_name,
                                            cols = sam_no_neg$accession)


obj_filt$data$sqrt_tax_abund <- calc_taxon_abund(obj_filt, "otu_sqrt",
                                            groups = NULL,#sam_no_neg$sample_name,
                                            cols = sam_no_neg$accession)
obj_filt$data$diff_table_tax_abund <- compare_groups(obj_filt,
                            dataset = "tax_abund",
                            cols = sam_no_neg$accession, # What columns of sample data to use
                            groups = sam_no_neg$sample_type)
obj_filt$data$diff_table_sqrt <- compare_groups(obj_filt,
                                                 dataset = "sqrt_tax_abund",
                                                 cols = sam_no_neg$accession, # What columns of sample data to use
                                                 groups = sam_no_neg$sample_type)

set.seed(3141)

heat_tree_matrix(obj_filt,
                 data = "diff_table_tax_abund",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-6,6),
                 edge_color_interval = c(-6,6),
                 node_size_axis_label = "XXXXXX transformed per-taxon counts",
                 node_color_axis_label = "Log2 ratio median diff.",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree.pdf") # Saves the plot as a pdf file


res<-full_join(obj_filt$data$tax_data,obj_filt$data$diff_table_tax_abund)

range(res$mean_diff)
range(res$median_diff)
range(res$log2_median_ratio)
favstats(res$log2_median_ratio)

write.table(res,"metacoder_supplementary_table.tsv",sep = "\t",row.names = F,quote = F)
node_
