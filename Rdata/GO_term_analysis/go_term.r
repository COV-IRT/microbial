library(phyloseq)
library(microbiome)
library(ggplot2)
library(microbiomeutilities)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggpubr)
library("pheatmap")
library(ggsci)
getwd()
count_tab<-read.table("Combined_BALF_GO_Terms_counts.txt", header = T, sep = "\t", row.names = 1)
tax_tab<-read.table("Combined_BALF_GO_Terms_tax.txt", header = T, sep = "\t", row.names = 1)
sam_tab<-read.table("Combined_BALF_GO_Terms_metadata.txt", header = T, sep = "\t", row.names = 1)
count<-otu_table(count_tab,taxa_are_rows = T)
tax<-tax_table(as.matrix(tax_tab),errorIfNULL = T)
sam<-sample_data(sam_tab)
pseq<-phyloseq(count,tax,sam)
summarize_phyloseq(pseq)

##########################################
mypal2<- pal_startrek("uniform", alpha = 1)(7)
mypal <- pal_aaas("default", alpha = 1)(7)
mypal2
mypal3<-rbind(mypal, mypal2)
#########################################

# Taxa with positive sum across samples
pseq_prune <- prune_taxa(taxa_sums(pseq) > 1, pseq)

min(sample_sums(pseq))
pseq_prune
library(ggsci)
library(RColorBrewer)
p <- plot_frequencies(sample_data(pseq), "publication","sample_type")+scale_fill_manual(values = mypal3)
print(p)

pseq_prune <- prune_taxa(taxa_sums(pseq) > 1, pseq)

rank_names(pseq_prune)
pseq_prune<-tax_glom(physeq = pseq_prune, taxrank = "name")
taxa_names(pseq_prune)<-get_taxa_unique(pseq_prune,taxonomic.rank = "name" )

dat <- abundances(pseq_prune)
count <- as.matrix(t(dat))
fit <- lapply(1:6, dmn, count = count, verbose=TRUE)

#Check the model fit with different number of mixture componenets using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
plot(aic)
plot(bic)
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)
best <- fit[[which.min(lplc)]]
best
best <-fit[[4]]
heatmapdmn(count, fit[[1]], best,ntaxa = 30,
           transform = sqrt, lblwidth = 0.2 * nrow(count))
mixturewt(best)
ass <- apply(mixture(best), 1, which.max)
ass
write.table(ass,"GO_TERMS_DMM.tsv",sep="")
write.table(fitted(best),"GO_term_DMMs.tsv", sep="\t")
physeq<-pseq_prune
meta<-data.frame(sample_data(physeq))
meta$dmm<-ass
for (k in seq(ncol(fitted(best)))) 
{
  d <- melt(fitted(best))
  colnames(d) <- c("GO", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange GOs by assignment strength
    arrange(value) %>%
    mutate(GO = factor(GO, levels = unique(GO))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = GO, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers in  : GO  Terms cluster type ", k, sep = ""))
  #paste(p,k, sep = "")<-p
  #print(k)
  print(p)
  #print(paste(p,k, sep=""))
}

rownames(ass)
colnames(ass)
sample_data(pseq_prune)$dmn<-ass


tmp<-psmelt(pseq_prune)
tmp<-as_tibble(tmp)
d2<-tmp%>%
  select(name,Abundance, dmn)%>%
  group_by(name, dmn) %>%
  summarise(avg = mean(Abundance)) %>%
  arrange(desc(avg))
d2<-d2%>%group_by(name, dmn)
d2
d3<-tidyr::spread(d2,dmn, avg)
d3[2:5]
d3$tot<-rowSums(d3[2:5])
d3
d3<-d3%>%arrange(desc(tot))
d3$tot<-NULL
d3<-d3[1:25,]
d3
d3<-d3%>%gather(data = d3,avg,2:5)
colnames(d3)<-c("name", "dmn","avg")
d3
d4<-data.frame(d3)
d4

ggballoonplot(d4, y ="name",x = "dmn", size = "avg", fill = "avg")+scale_fill_viridis_c()


count<-abundances(pseq_prune)
count
rowMeans(count,na.rm = T)

select <- order(rowMeans(count),
                decreasing=TRUE)[1:20]
select
#sample_data(pseq_prune)<-sam
dim(sam)
dim(sample_data(pseq_prune))
df<-as.data.frame(sam)
df<-as_tibble(df)
df<-df%>%select(dmn)#,dmn,body_site)
df <-sample_data(pseq_prune)$dmn

df<-as.character(df)
df<-as.data.frame(df)
df$Dataset<-sam$publication

df$sample_type<-sample_data(pseq_prune)$sample_type
select2<-(count)[select,]
colnames(select2) <- colnames(otu_table(pseq_prune))
colnames(select2)
row.names(df)
row.names(df) <- colnames(select2)
ass
mypal <- pal_aaas("default", alpha = 1)(10)
mypal
library("scales")
unique(df)
show_col(mypal)
df
mosaicplot(df)
# Specify colors
ann_colors = list(
  df=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF"),
  Dataset=c("Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="#631879FF",       
  "Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF",
  "Huang"="#000000","Ren"="#111111"),
  sample_type =c("COVID_19"="#3B4992FF",
                 "Community_acquired_pneumonia"="#EE0000FF", 
                 "Healthy"="#008B45FF",                     
                  "neg_control"="#000000",
                 "Smoker"="#631879FF",
                 "Asthma"= "#008280FF",                      
                  "Asthma_Smoker"="#BB0021FF",
                 "Asthma_Ex_smoker"="#5F559BFF",
                 "Obese"="#A20056FF",                       
                  "Obese_Asthma"="#808180FF",
                 "Obese_Smoker"="#1B1919FF",
                 "Obese_Asthma_Smoker"="#1111111",
                 "Unknown"="#222222"))
ann_colors
df_row<-as.data.frame(fitted(best))
colnames(df_row)L-c()
library(RColorBrewer)
library(viridis)
xx <- pheatmap(mat = select2,
               color = viridis(256),
               annotation_col=df, 
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               annotation_row = df_row)
xx
pheatmap(mat, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
         kmeans_k = NA, 
         breaks = NA, 
         border_color = "grey60",
         cellwidth = NA, 
         cellheight = NA, 
         scale = "none", 
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         clustering_callback = identity2, 
         cutree_rows = NA, 
         cutree_cols = NA,
         treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
                                 50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
                                                                   cluster_cols, 50, 0), 
         legend = TRUE,
         legend_breaks = NA,
         legend_labels = NA,
         annotation_row = NA, 
         annotation_col = NA,
         annotation = NA, 
         annotation_colors = NA, 
         annotation_legend = TRUE,
         annotation_names_row = TRUE, 
         annotation_names_col = TRUE,
         drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
         fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
         angle_col = c("270", "0", "45", "90", "315"), display_numbers = F,
         number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8
         * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,
         silent = FALSE, na_col = "#DDDDDD", ...)