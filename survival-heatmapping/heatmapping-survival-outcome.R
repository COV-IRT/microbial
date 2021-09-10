library(phyloseq)
library(microbiome)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)

setwd("github/microbial/survival-heatmapping/")
# save.image()
# load(".RData") # this comes from the original heatmapping.R file
# load("2_maaslin2.rda")
  ## code is modifying that initial code (initial code is commented out)
# compsitionally transform both datasets

load("data-for-outcomes-heatmap-example.RData") # this is what was saved at the end, after building off of the above
term_pseq_no_neg_gonames
# sample_data(term_pseq_no_neg_gonames)$sample_type
# term_pseq_no_neg_gonames_COVID19<-subset_samples(physeq = term_pseq_no_neg_gonames,sample_type=="COVID_19")
term_pseq_comp <- microbiome::transform(transform = "compositional", term_pseq_no_neg_gonames)
#term_pseq_comp <- microbiome::transform(transform = "compositional", term_pseq_no_neg_gonames_COVID19)
term_pseq_comp
GO_terms_to_plot <- scan("GO-terms-to-plot.txt", what="character")

# term_pseq_prune_comp2<-microbiome::transform(transform = "compositional",term_pseq_no_neg_gonames)
# term_pseq_prune_comp<-microbiome::transform(transform = "compositional",term_pseq_prune)
#
# # prune out the masslin2 GO Terms and fix the names for the big dataset
term_pseq_prune_comp <- prune_taxa(taxa = GO_terms_to_plot, term_pseq_comp)
tax <- data.frame(tax_table(term_pseq_prune_comp))
names<-paste(rownames(tax), tax$name, sep = "-")
taxa_names(term_pseq_prune_comp) <- names

#
# term_pseq_prune_comp2 <- prune_taxa(taxa = Terms$Term,x =term_pseq_prune_comp2)
# tax<-data.frame(tax_table(term_pseq_prune_comp2))
# names<-paste(rownames(tax),tax$name,sep="-")
# taxa_names(term_pseq_prune_comp2)<-names

# order the count tables by the rowMeans
#
# count<-abundances(term_pseq_prune_comp)
# select <- order(rowMeans(count),decreasing=TRUE)
# select2<-(count)[select,]
# tmp<-rownames(select2)
# select2<-as_tibble(select2, rownames=NA)
#
# countb<-abundances(term_pseq_prune_comp2)
# selectb <- order(rowMeans(countb),decreasing=TRUE)
# select2b<-(countb)[select,]
# tmpb<-rownames(select2b)
# select2b<-as_tibble(select2b, rownames=NA)
# head(colnames(select2))
#
#
# # make an annotation dataset for the colum headers
# sam<-data.frame(sample_data(term_pseq_prune_comp))
# df<-as.data.frame(sample_data(term_pseq_prune_comp))
# all.equal(row.names(df), colnames(select2))
#
# df<-as_tibble(df, rownames = NA)
# df<-df%>%select(publication, dmn,case,outcome)
# df<-as.data.frame(df)
# row.names(df) <- colnames(select2)
# colnames(df)<-c("Publication","dmm_cluster", "Case","Outcome")
# df$dmm_cluster<-as.character(df$dmm_cluster)
#
#   ## naming this one as samb (even if some of this is the same, it gets confusing and more prone to mistake if using the same variable names but different code)
# samb<-data.frame(sample_data(term_pseq_prune_comp2))
# dfb<-as.data.frame(sample_data(term_pseq_prune_comp2))
# dfb<-as_tibble(dfb, rownames = NA)
# dfb<-dfb%>%select(publication,dmn, case,outcome)
# dfb<-as.data.frame(dfb)
# row.names(dfb) <- colnames(select2b)
# colnames(dfb)<-c("Publication","dmm_cluster", "Case","Outcome")
# dfb$dmm_cluster<-as.character(dfb$dmm_cluster)
#
# # make the annotation rows
#
# df_row<-as.data.frame(fitted(best))
#   ## don't know if these new column names are ideal yet, or if just adding something later would help. but one other way to name these in the legend
# colnames(df_row)<-c("Contr. to DMM 1","Contr. to DMM 2","Contr. to DMM 3")
# # make the associated color pallete for the column and row headers
# brewer.pal(9, "Paired")
#
# ann_colors = list(
#     Publication=c("Xiong"="#008B45FF", "Shen"="#3B4992FF","Chen"='#631879FF',"Wu"="#008280FF","Zhou"='#BB0021FF',"Ren"='#FF7F00'),
#     dmm_cluster=c("1"="forestgreen","2"="darkorange1","3"="firebrick"),
#     Case=c("COVID19"="firebrick","Community_acquired_pneumonia"="darkorange1", "Control_Healthy"="forestgreen"),
#     Outcome=c("Deceased"="black","Stabilized"="goldenrod1","Recovered"="forestgreen"))
#
# dim(df_row)
#
# subset_row_annots <- df_row[rownames(df_row) %in% rownames(select2b), ] # hmm, why are none of these in here
#
# head(rownames(select2b))
# # what is "GO:0003824-catalytic activity" labeled as in df_row?
#
# rownames(df_row)[rownames(df_row) == "GO:0003824-catalytic activity"]
#   # nothing
#
# rownames(df_row)[grep("GO:0003824", rownames(df_row))]
#   # "GO:0003824-catalytic activity-catalytic activity", hmm
# wanted_ids <- vector()
# for ( id in rownames(select2b) ) {
#     wanted_ids <- c(wanted_ids, rownames(df_row)[grep(id, rownames(df_row))])
# }
#
# subset_row_annots <- df_row[rownames(df_row) %in% wanted_ids, ]
# summary(subset_row_annots)
#
#   # df_row was made from this fitted(best)
# head(fitted(best))
#     # which does have all their names twice, e.g.:
#     # "GO:0000002-mitochondrial genome maintenance-mitochondrial genome maintenance"
#
#   # not an easy fix here because it'll require going back to how "best" was made, so dropping the dmm weighting column
#   # see what you think about normalizing them too, like making them all percents of their totals, e.g.:
#
# mod_row_annots_df <- df_row %>% apply(2, function(x) ( x / sum(x)) * 100)
#     # just cause it don't think it's very interpretable having them being different colors and different scales,
#     # as it is, I don't we can look at one GO term and see which DMM it contributed to most, which i think is the purpose of having this there?
#
# xx <- pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                annotation_col=df,
#                cutree_col = 14,
#                cutree_row=3,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors,
#                clustering_distance_rows = "euclidean",
#                clustering_distance_cols = "euclidean")
#
# xx <- pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                annotation_col=df,
#                cutree_col = 14,
#                cutree_row=3,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors,
#                clustering_distance_rows = "euclidean",
#                clustering_distance_cols = "euclidean",
#                annotation_row=df_row)
#
# xx <- pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                annotation_col=df,
#                cutree_col = 14,
#                cutree_row=3,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors,
#                clustering_distance_rows = "euclidean",
#                clustering_distance_cols = "euclidean",
#                annotation_row=subset_row_annots)
#
#


##### making our own clustering #####
  # trying with all, VST -> euclidean dist -> ward.D2 clustering just because
deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(term_pseq_no_neg_gonames)), colData = data.frame(sample_data(term_pseq_no_neg_gonames)), design = ~1)

#deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(term_pseq_no_neg_gonames_COVID19)), colData = data.frame(sample_data(term_pseq_no_neg_gonames_COVID19)), design = ~1)
vst_counts <- varianceStabilizingTransformation(deseq_counts)
vst_counts_df <- assay(vst_counts)

sample_euc_dist <- dist(t(vst_counts_df))
# not making a GO/row one because we won't have all of the terms in the figure, and we want to order them our own way anyway
sample_hclust <- hclust(sample_euc_dist, method = "ward.D2")
plot(sample_hclust)

# pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                cluster_cols = sample_hclust,
#                annotation_col=df,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors)
#
#   # blocking out the 4 main sample Hclusters
# pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                cluster_cols = sample_hclust,
#                cutree_cols = 4,
#                annotation_col=df,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors)
#
#   # no tree on GO/row
# pheatmap(mat = select2b,
#                color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(20),
#                cluster_cols = sample_hclust,
#                cluster_rows = FALSE,
#                cutree_cols = 4,
#                annotation_col=df,
#                scale="row",
#                border_color="NA",
#                annotation_colors = ann_colors)


##############################
# grouped sig GO terms in "maaslin2-outcome-survival-heatmap.xlsx", table we're reading in is from the tab "maaslin2-outcome-heatmap-info"
##############################

GO_annots_heatmap_info <- read.table("maaslin2-outcome-survival-heatmap-plotting-info.tsv", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
# this object holds the new table I put together for grouping the GO terms in the heatmap
str(GO_annots_heatmap_info)

percentages_df <- data.frame(abundances(term_pseq_prune_comp))
percentages_df
# reading in starting GO info tab
# GO_annots_heatmap_info$`Depth 1 Parent(s)` <- factor(GO_annots_heatmap_info$`Depth 1 Parent(s)`, levels = unique(GO_annots_heatmap_info$`Depth 1 Parent(s)`))
# subsetting down to those in the GO_annots_heatmap_info tab
# select2b_df <- data.frame(select2b)
# sub_GO_for_heatmap_vis_df <- select2b_df[row.names(select2b_df) %in% row.names(GO_annots_heatmap_info), ]

# ordering plot
# sub_GO_for_heatmap_vis_df <- sub_GO_for_heatmap_vis_df[row.names(GO_annots_heatmap_info), ]
GO_for_heatmap <- percentages_df[row.names(GO_annots_heatmap_info), ]

# gaps vector
gaps_vec <- c(3, 12, 16)

# making color vector for Depth 1 parents col

library(cartography)
GO_annots_heatmap_info$`Depth 1 Parent(s)`<-as_factor(GO_annots_heatmap_info$`Depth 1 Parent(s)`)
length(levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`))

parent_cols <- carto.pal("multi.pal", n1 = length(levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`)))
names(parent_cols) <- levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`)

# colors for GO NameSpace col
namespace_cols <- c("Biological Process" = "#3B4992FF", "Molecular Function" = "#008B45FF")

# new annotation colors obj.
annot_cols <- c(ann_colors, list("Depth 1 Parent(s)" = parent_cols, "NameSpace" = namespace_cols))

## adding a grey color and "NA" label for Outcome for those without any, so it's not empty on the plot
sample_annots_df
# 
# sample_annots_df2<-sample_annots_df%>%filter(Case=="COVID19")
# GO_for_heatmap2<-GO_for_heatmap%>%select(rownames(sample_annots_df2))
# GO_for_heatmap2
#lets chack the stabilized and recoved to Survived and fix the colors
sample_annots_df$Outcome<-gsub("Recovered","Survived",sample_annots_df$Outcome)
sample_annots_df$Outcome<-gsub("Stabilized","Survived",sample_annots_df$Outcome)

sample_annots_df$Case<-gsub("Control_Healthy","Uninfected",sample_annots_df$Case)
#sample_annots_df$Outcome
# sample_annots_df$Outcome[is.na(sample_annots_df$Outcome)] <- "NA"
annot_cols$Outcome <- c("Deceased" = "black","Survived" = "goldenrod1","NA" = "grey90")
#annot_cols$Outcome <- c(annot_cols$Outcome, "NA" = "grey90")
annot_cols$Case

annot_cols$Case <- c("COVID19" = "firebrick","Community_acquired_pneumonia" = "darkorange1","Uninfected" = "forestgreen")

sample_hclust

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
               cellwidth = 15, cellheight = 15, fontsize = 10,
               filename = "GO-outcome-survival-heatmap-only-those-with-outcomes.pdf", width = 22, height = 9, main = "Significantly different GO Terms")

dev.off()
getwd()
save.image("data-for-outcomes-heatmap-example.RData")

