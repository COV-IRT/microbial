#load the necessady libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
#set the working directory
setwd("c:/github/microbial/GO_Terms_scripts/scripts/Lee/")
#load the necessary tables and phyloseq images
#phyloseq image
term_pseq_no_neg_gonames<-readRDS("term_pseq_no_neg_gonames.RDS")
#names of the GO terms to plot
GO_terms_to_plot <- scan("c:/github/microbial/survival-heatmapping/GO-terms-to-plot.txt", what="character")
#heatmap annotation info
GO_annots_heatmap_info <- read.table("../../../survival-heatmapping/maaslin2-outcome-survival-heatmap-plotting-info.tsv", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
sample_annots_df<-read.table(file = "sample_annots_df.tsv",header = T,sep = "\t",row.names = 1)

#conver the read counts to relative abundances
term_pseq_comp<-microbiome::transform(term_pseq_no_neg_gonames,"compositional")

# # prune out the masslin2 GO Terms and fix the names for the big dataset
term_pseq_prune_comp <- prune_taxa(taxa = GO_terms_to_plot, term_pseq_comp)
tax <- data.frame(tax_table(term_pseq_prune_comp))
names<-paste(rownames(tax), tax$name, sep = "-")
taxa_names(term_pseq_prune_comp) <- names



##### making our own clustering #####
# trying with all, VST -> euclidean dist -> ward.D2 clustering just because
deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(term_pseq_no_neg_gonames)), colData = data.frame(sample_data(term_pseq_no_neg_gonames)), design = ~1)

#deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(term_pseq_no_neg_gonames_COVID19)), colData = data.frame(sample_data(term_pseq_no_neg_gonames_COVID19)), design = ~1)
vst_counts <- varianceStabilizingTransformation(deseq_counts)
vst_counts_df <- assay(vst_counts)

sample_euc_dist <- dist(t(vst_counts_df))
# not making a GO/row one because we won't have all of the terms in the figure, and we want to order them our own way anyway
sample_hclust <- hclust(sample_euc_dist, method = "ward.D2")

##############################
# grouped sig GO terms in "maaslin2-outcome-survival-heatmap.xlsx", table we're reading in is from the tab "maaslin2-outcome-heatmap-info"
##############################
# this object holds the new table I put together for grouping the GO terms in the heatmap
str(GO_annots_heatmap_info)

percentages_df <- data.frame(abundances(term_pseq_prune_comp))
# ordering plot
GO_for_heatmap <- percentages_df[row.names(GO_annots_heatmap_info), ]

# gaps vector
gaps_vec <- c(3, 12, 16)
length(levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`))


####color blind fiendly edits
# make the associated color pallete for the column and row headers
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot_cols = list(Publication=c("Xiong"="#E69F00","Shen"="#56B4E9","Chen"= "#009E73",
                                "Wu"="#F0E442","Zhou"="#0072B2","Ren"= "#D55E00"),
                  dmm_cluster=c("1"="#0072B2","2"= "#D55E00","3"="#CC79A7"),
                  Case=c("COVID19"= "#D55E00","Community_acquired_pneumonia"="#F0E442","Uninfected"="#009E73"),
                  Outcome=c("Deceased"="#000000","Survived"="#F0E442","NA"="#999999"),
                  Namespace= c("Biological Process" = "#071f63", "Molecular Function" = "#b47b00"),
                  "Depth 1 Parent(s)" =c("GO:0008152 | metabolic process"="#b28c3b",
                                         "GO:0008152 metabolic process | GO:0009987  cellular process"="#8364b8",
                                         "GO:0003824 | catalytic activity"="#57af6c","GO:0005488 | binding"="#b74665"))


#write.table(sample_annots_df,"sample_annots_df.tsv",sep = "\t",row.names = F)
#
#lets chack the stabilized and recoved to Survived and fix the colors
sample_annots_df$Outcome<-gsub("Recovered","Survived",sample_annots_df$Outcome)
sample_annots_df$Outcome<-gsub("Stabilized","Survived",sample_annots_df$Outcome)
sample_annots_df$Case<-gsub("Control_Healthy","Uninfected",sample_annots_df$Case)

#edited for colorblind friendly


pheatmap(mat = GO_for_heatmap,
               color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
               cluster_cols = sample_hclust,
               cluster_rows = FALSE,
               cutree_cols = 4,
               annotation_col = sample_annots_df,
               annotation_row = GO_annots_heatmap_info,
               scale="row",
               border_color="black",
               annotation_colors = annot_cols,
               angle_col = 315,
               fontsize_row = 8, fontsize_col = 7, gaps_row = gaps_vec,
               cellwidth = 10, cellheight = 10, fontsize = 9)

# writing one out
pheatmap(mat = GO_for_heatmap,
               color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
               cluster_cols = sample_hclust,
               cluster_rows = FALSE,
               cutree_cols = 4,
               annotation_col = sample_annots_df,
               annotation_row = GO_annots_heatmap_info,
               scale="row",
               border_color="black",
               annotation_colors = annot_cols,
               angle_col = 315,
               fontsize_row = 8, fontsize_col = 7, gaps_row = gaps_vec,
               cellwidth = 10, cellheight = 10, fontsize = 9,
               filename = "GO-outcome-survival-heatmap.pdf", width = 22, height = 9, main = "Significantly different GO Terms")


dev.off()

###### making one with just those that have outcomes ######
samples_with_outcomes_sample_annots_df <- sample_annots_df %>% filter(Outcome != "NA")
#samples_with_outcomes_sample_annots_df <- sample_annots_df %>% filter(Case == "COVID19")
dim(samples_with_outcomes_sample_annots_df)
samples_with_outcomes_sample_annots_df <- samples_with_outcomes_sample_annots_df %>% select(- Case)

samples_with_outcomes_term_pseq_no_neg_gonames <- subset_samples(term_pseq_no_neg_gonames, outcome %in% samples_with_outcomes_sample_annots_df$Outcome %>% unique())

samples_with_outcomes_term_pseq_no_neg_gonames <- subset_samples(term_pseq_no_neg_gonames, !is.na(outcome))
samples_with_outcomes_term_pseq_no_neg_gonames
#samples_with_outcomes_term_pseq_no_neg_gonames <- subset_samples(term_pseq_no_neg_gonames, case == "COVID19")

samples_with_outcomes_term_pseq_no_neg_gonames 
# new clustering (needed because we removed samples)
samples_with_outcomes_deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(samples_with_outcomes_term_pseq_no_neg_gonames)), colData = data.frame(sample_data(samples_with_outcomes_term_pseq_no_neg_gonames)), design = ~ publication)
samples_with_outcomes_vst_counts <- varianceStabilizingTransformation(samples_with_outcomes_deseq_counts)
samples_with_outcomes_vst_counts_df <- assay(samples_with_outcomes_vst_counts)

samples_with_outcomes_sample_euc_dist <- dist(t(samples_with_outcomes_vst_counts_df))
# not making a GO/row one because we won't have all of the terms in the figure, and we want to order them our own way anyway
samples_with_outcomes_sample_hclust <- hclust(samples_with_outcomes_sample_euc_dist, method = "ward.D2")


samples_with_outcomes_percentages_df <- percentages_df %>% select(row.names(samples_with_outcomes_sample_annots_df))

samples_with_outcomes_GO_for_heatmap <- samples_with_outcomes_percentages_df[row.names(GO_annots_heatmap_info), ]


pheatmap(mat = samples_with_outcomes_GO_for_heatmap,
               color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
               cluster_cols = samples_with_outcomes_sample_hclust,
               cluster_rows = FALSE,
               cutree_cols = 2,
               annotation_col = samples_with_outcomes_sample_annots_df,
               annotation_row = GO_annots_heatmap_info,
               scale="row",
               border_color="black",
               annotation_colors = annot_cols,
               angle_col = 315,
               fontsize_row = 10, fontsize_col = 9, gaps_row = gaps_vec,
               cellwidth = 15, cellheight = 15, fontsize = 10)

# writing one out
pheatmap(mat = samples_with_outcomes_GO_for_heatmap,
               #color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
               cluster_cols = samples_with_outcomes_sample_hclust,
               cluster_rows = FALSE,
               cutree_cols = 2,
               annotation_col = samples_with_outcomes_sample_annots_df,
               annotation_row = GO_annots_heatmap_info,
               scale="row",
               border_color="black",
               annotation_colors = annot_cols,
               angle_col = 315,
               fontsize_row = 10, fontsize_col = 9, gaps_row = gaps_vec,
               cellwidth = 15, cellheight = 15, fontsize = 10,
               filename = "Figure3_GO-outcome-survival-heatmap-only-those-with-outcomes_colorblind_friendly.pdf", 
         width = 22, height = 9, main = "Significantly different GO Terms")

dev.off()
getwd()
save.image("data-for-outcomes-heatmap-example.RData")

