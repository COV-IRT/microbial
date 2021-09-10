
##########
## This takes off from "heatmapping.R" we shared a while back in slack channel
##########

setwd("~/github/microbial/")
# this loads a subset of the data needed to do this part
load("sub-data-for-heatmap-example.RData")

# this object holds the new table I put together for grouping the GO terms in the heatmap
head(GO_annots_heatmap_info)

# subsetting former GO heatmap down to those in the GO_annots_heatmap_info tab (i removed those that are level 1s, see tab 3 of newest "sig-go-info.xlsx")
select2b_df <- data.frame(select2b)
sub_GO_for_heatmap_vis_df <- select2b_df[row.names(select2b_df) %in% row.names(GO_annots_heatmap_info), ]

# ordering plot (factor business isn't needed for this)
sub_GO_for_heatmap_vis_df <- sub_GO_for_heatmap_vis_df[row.names(GO_annots_heatmap_info), ]

# gaps vector (horizontal spaces we want in the plot)
gaps_vec <- c(4, 5, 13, 18, 20)

# making color vector for Depth 1 parents col
library(cartography) # vector is pasted below if you need the colors and don't want to install cartography
parent_cols <- carto.pal("multi.pal", n1 = length(levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`)))
names(parent_cols) <- levels(GO_annots_heatmap_info$`Depth 1 Parent(s)`)
parent_cols

#    GO:0003824 | catalytic activity                                        GO:0005488 | binding
#                          "#4e3f41"                                                   "#76c663"
#     GO:0008152 | metabolic process GO:0008152 metabolic process | GO:0009987  cellular process
#                          "#9451b2"                                                   "#bda856"
# GO:0065007 | biological regulation     GO:0044419 | interspecies interaction between organisms
#                          "#92afbd"                                                   "#c4544e"

# colors for GO NameSpace col
namespace_cols <- c("Biological Process" = "#3B4992FF", "Molecular Function" = "#008B45FF")

# new annotation colors object
annot_cols <- c(ann_colors, list("Depth 1 Parent(s)" = parent_cols, "NameSpace" = namespace_cols))

## adding a grey color and "NA" label for Outcome for those without any, so it's not empty on the plot
sample_annots_df <- df
sample_annots_df$Outcome[is.na(sample_annots_df$Outcome)] <- "NA"

annot_cols$Outcome <- c(annot_cols$Outcome, "NA" = "grey90")

library(pheatmap)
library(RColorBrewer)
pheatmap(mat = sub_GO_for_heatmap_vis_df,
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
pheatmap(mat = sub_GO_for_heatmap_vis_df,
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
         filename = "GO-heatmap.pdf", width = 22, height = 9, main = "Significantly different GO Terms")


save.image("sub-data-for-heatmap-example.RData")
