# Figure 3 - Outcome heatmapping plot


#load the necessady libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
#set the working directory
setwd("K:/github/microbial/Figures/02092022/Figure3/")
#load the necessary tables and phyloseq images
#phyloseq image
term_pseq_no_neg_gonames<-readRDS("term_pseq_no_neg_gonames.RDS")
#names of the GO terms to plot
GO_terms_to_plot <- scan("GO-terms-to-plot.txt", what="character")
#heatmap annotation info
GO_annots_heatmap_info <- read.table("maaslin2-outcome-survival-heatmap-plotting-info.tsv", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
sample_annots_df<-read.table(file = "sample_annots_df.tsv",header = T,sep = "\t",row.names = 1)
#write.table(sample_annots_df,file = "sample_annots_df.tsv",sep="\t",row.names = T)

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



sample_euc_dist <- dist(t(data.frame(otu_table(term_pseq_no_neg_gonames))))
sample_euc_dist

sample_euc_dist <- dist(t(vst_counts_df))

# not making a GO/row one because we won't have all of the terms in the figure, and we want to order them our own way anyway
sample_hclust <- hclust(sample_euc_dist, method = "ward.D2")


# Plot the obtained dendrogram
plot(sample_hclust, cex = 0.6, hang = -1)
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



####color blind fiendly edits
# make the associated color pallete for the column and row headers
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

annot_cols = list(Publication=c("Xiong"="#E69F00","Shen"="#56B4E9","Chen"= "#009E73",
                                "Wu"="#F0E442","Zhou"="#0072B2","Ren"= "#D55E00"),
                  dmm_cluster=c("1"="#0072B2","2"= "#D55E00","3"="#CC79A7"),
                  Case=c("COVID19"= "#D55E00","Community_acquired_pneumonia"="#F0E442","Uninfected"="#009E73"),
                  Outcome=c("Deceased"="#000000","Survived"="#F0E442","NA"="#999999"),
                  Namespace= c("Biological Process" = "#071f63", "Molecular Function" = "#b47b00"),
                  "Depth 1 Parent(s)" =c("GO:0008152 | metabolic process"="#b28c3b",
                                         "GO:0008152 metabolic process | GO:0009987  cellular process"="#8364b8",
                                         "GO:0003824 | catalytic activity"="#57af6c","GO:0005488 | binding"="#b74665"))




# pheatmap(mat = GO_for_heatmap,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                cluster_cols = sample_hclust,
#                cluster_rows = FALSE,
#                cutree_cols = 4,
#                annotation_col = sample_annots_df,
#                annotation_row = GO_annots_heatmap_info,
#                scale="row",
#                border_color="black",
#                annotation_colors = annot_cols,
#                angle_col = 315,
#                fontsize_row = 8, fontsize_col = 7, gaps_row = gaps_vec,
#                cellwidth = 10, cellheight = 10, fontsize = 9)
# 
# # writing one out
# pheatmap(mat = GO_for_heatmap,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                cluster_cols = sample_hclust,
#                cluster_rows = FALSE,
#                cutree_cols = 4,
#                annotation_col = sample_annots_df,
#                annotation_row = GO_annots_heatmap_info,
#                scale="row",
#                border_color="black",
#                annotation_colors = annot_cols,
#                angle_col = 315,
#                fontsize_row = 8, fontsize_col = 7, gaps_row = gaps_vec,
#                cellwidth = 10, cellheight = 10, fontsize = 9,
#                filename = "GO-outcome-survival-heatmap.pdf", width = 22, height = 9, main = "Significantly different GO Terms")
# 
# 
# dev.off()

###### making one with just those that have outcomes ######
samples_with_outcomes_sample_annots_df <- sample_annots_df %>% filter(Outcome != "NA")
samples_with_outcomes_sample_annots_df <- samples_with_outcomes_sample_annots_df %>% select(- Case)
samples_with_outcomes_term_pseq_no_neg_gonames <- subset_samples(term_pseq_no_neg_gonames, outcome %in% samples_with_outcomes_sample_annots_df$Outcome %>% unique())
samples_with_outcomes_term_pseq_no_neg_gonames <- subset_samples(term_pseq_no_neg_gonames, !is.na(outcome))
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


round(range(samples_with_outcomes_GO_for_heatmap),2)
#new annotations list with only the publications used in the analysis
# annot_cols = list(Publication=c("Shen"="#56B4E9","Chen"= "#009E73","Zhou"="#0072B2"),
#                   dmm_cluster=c("1"="#0072B2","2"= "#D55E00","3"="#CC79A7"),
#                   Case=c("COVID19"= "#D55E00","Community_acquired_pneumonia"="#F0E442","Uninfected"="#009E73"),
#                   Outcome=c("Deceased"="#000000","Survived"="#F0E442","NA"="#999999"),
#                   Namespace= c("Biological Process" = "#071f63", "Molecular Function" = "#b47b00"),
#                   "Depth 1 Parent(s)" =c("GO:0008152 | metabolic process"="#b28c3b",
#                                          "GO:0008152 metabolic process | GO:0009987  cellular process"="#8364b8",
#                                          "GO:0003824 | catalytic activity"="#57af6c","GO:0005488 | binding"="#b74665"))

annot_cols = list(Publication=c("Shen"="#56B4E9","Chen"= "#009E73","Zhou"= "#D55E00"),
                  dmm_cluster=c("1"="#0072B2","2"= "#D55E00","3"="#CC79A7"),
                  Case=c("COVID19"= "#D55E00","Community_acquired_pneumonia"="#F0E442","Uninfected"="#009E73"),
                  Outcome=c("Deceased"="#000000","Survived"="#F0E442","NA"="#999999"),
                  NameSpace= c("Biological Process" = "#0C7BDC", "Molecular Function" = "#FFC20A"),
                  "Depth 1 Parent(s)" =c("GO:0008152 | metabolic process"="#E69F00",
                                         "GO:0008152 metabolic process | GO:0009987  cellular process"="#56B4E9",
                                         "GO:0003824 | catalytic activity"="#009E73",
                                         "GO:0005488 | binding"="#F0E442"))
meths<-c('ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' , 'centroid')
library(vegan)



# new clustering (needed because we removed samples)

library(tidyverse)
sam<-data.frame(sample_data(samples_with_outcomes_term_pseq_no_neg_gonames))
samples_with_outcomes_GO_for_heatmap_in<-samples_with_outcomes_GO_for_heatmap*10000
range(samples_with_outcomes_GO_for_heatmap_in)
class(samples_with_outcomes_GO_for_heatmap_in)
samples_with_outcomes_sample_annots_df
samples_with_outcomes_deseq_counts <- DESeqDataSetFromMatrix(samples_with_outcomes_GO_for_heatmap_in[1],
                                                             colData = samples_with_outcomes_sample_annots_df, 
                                                             design = ~ Outcome)
samples_with_outcomes_vst_counts <- varianceStabilizingTransformation(samples_with_outcomes_deseq_counts)
samples_with_outcomes_vst_counts_df <- assay(samples_with_outcomes_vst_counts)
samples_with_outcomes_sample_euc_dist <- dist(t(samples_with_outcomes_vst_counts_df))
# not making a GO/row one because we won't have all of the terms in the figure, and we want to order them our own way anyway
samples_with_outcomes_sample_hclust <- hclust(samples_with_outcomes_sample_euc_dist, method = "ward.D2")
samples_with_outcomes_percentages_df <- percentages_df %>% select(row.names(samples_with_outcomes_sample_annots_df))
samples_with_outcomes_GO_for_heatmap <- samples_with_outcomes_percentages_df[row.names(GO_annots_heatmap_info), ]

c<-hclust(vegdist(x = samples_with_outcomes_GO_for_heatmap,method = "manhattan"))


pheatmap(mat = samples_with_outcomes_GO_for_heatmap,
               color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#               cluster_cols =samples_with_outcomes_sample_hclust,
         cluster_cols =c,
         #clustering_method = meths[5],
               cluster_rows = FALSE,
               cutree_cols = 2,
               annotation_col = samples_with_outcomes_sample_annots_df,
               annotation_row = GO_annots_heatmap_info,
               scale="row",
         #main = paste0(meths[5]),
               border_color="black",
               annotation_colors = annot_cols,
               angle_col = 315,
               fontsize_row = 10, fontsize_col = 9, gaps_row = gaps_vec,
               cellwidth = 15, cellheight = 15, fontsize = 10)

# writing one out
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
               fontsize_row = 20, fontsize_col = 18, gaps_row = gaps_vec,
               cellwidth = 30, cellheight = 30, fontsize = 20,
              # filename = "Figure3_GO-outcome-survival-heatmap-only-those-with-outcomes_colorblind_friendly.pdf", 
         width = 32, height = 14, main = "Significantly different GO Terms")

dev.off()
samples_with_outcomes_sample_hclust$dist.method

save.image("Figure3_survival_heatmap_02092022.RData")

