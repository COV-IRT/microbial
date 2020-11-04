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
# getwd()
# count_tab<-read.table("Combined_BALF_GO_Terms_counts.txt", header = T, sep = "\t", row.names = 1)
# tax_tab<-read.table("Combined_BALF_GO_Terms_tax.txt", header = T, sep = "\t", row.names = 1)
# sam_tab<-read.table("Combined_BALF_GO_Terms_metadata.txt", header = T, sep = "\t", row.names = 1)
# count<-otu_table(count_tab,taxa_are_rows = T)
# tax<-tax_table(as.matrix(tax_tab),errorIfNULL = T)
# sam<-sample_data(sam_tab)
pseq<-bio_bac_physeq
pseq<-subset_samples(pseq, sample_type!="neg_control")
pseq<-subset_samples(pseq, sample_type!="Unknown")

summarize_phyloseq(pseq)
meta(pseq)$case
##########################################
mypal2<- pal_startrek("uniform", alpha = 1)(7)
mypal <- pal_aaas("default", alpha = 1)(7)
mypal2
mypal3<-rbind(mypal, mypal2)
mypal3
#########################################

# Taxa with positive sum across samples
pseq_prune <- prune_taxa(taxa_sums(pseq) > 1, pseq)
pseq_prune
min(sample_sums(pseq))
pseq_prune

mypal3<-c("#3B4992FF", "#BB0021FF", "#5F559BFF", "#CC0C00FF", "#EE0000FF","#008B45FF",  "#631879FF", "#008280FF", "#5C88DAFF", "#84BD00FF", "#FFCD00FF", "#7C878EFF", "#00B5E2FF", "#00AF66FF")
library(ggsci)
library(RColorBrewer)
p <- plot_frequencies(sample_data(pseq), "publication","sample_type")+scale_fill_manual(values = mypal3)
print(p)
q <- plot_frequencies(sample_data(pseq), "publication","case")+scale_fill_manual(values =c( "#008B45FF","#3B4992FF" ,"#EE0000FF"))#green-blue-red 
print(q)
pseq_prune <- prune_taxa(taxa_sums(pseq) > 10000, pseq)
pseq_prune <- prune_samples(sample_sums(pseq_prune) > 10, pseq_prune)
pseq_prune

pseq_prune<-tax_glom(physeq = pseq_prune, taxrank = "name")
taxa_names(pseq_prune)<-get_taxa_unique(pseq_prune,taxonomic.rank = "name" )
pseq_prune
dat <- abundances(pseq_prune)
count <- as.matrix(t(dat))
fit <- mclapply(1:6, dmn, count = count, verbose=TRUE)
fit
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
#best <-fit[[4]]
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
#tmp$OTU



tmp<-as_tibble(tmp)
tmp

d2<-tmp %>%
  select(case,name,Abundance, dmn)%>%
  group_by(name,case, dmn) %>%
  summarise(avg = mean(Abundance)) %>%
  arrange(desc(avg))

d2

d2
d2<-d2%>%group_by(name, case,dmn)
d2
d3<-tidyr::spread(d2,dmn, avg)
d3
d3[3:8]
d3$tot<-rowSums(d3[3:8])
d3
d3<-d3%>%arrange(desc(tot))
d3$tot<-NULL
d3<-d3[1:25,]
d3
d3<-d3%>%gather(data = d3,avg,3:8)
d3
colnames(d3)<-c("name","case", "dmn","avg")
d3
d4<-data.frame(d3)
d4

ggballoonplot(d4, y ="name",x = "dmn", size = "avg", fill = "avg")+scale_fill_viridis_c()
library(mosaic)

meta(pseq_prune)
head(tmp$OTU)
head(tmp$Abundance)
head(tmp$dmn)





my_tbl<-tally(~case+ dmn+publication,data =meta(pseq_prune),format = 'count')
my_tbl
a<-as_tibble(my_tbl)
a
ggballoonplot(data = a, y ="case",x = "dmn", facet.by = "publication",size = "n", fill = "n")+scale_fill_viridis_c(option = "C")
my_tbl<-tally(publication ~ dmn,data =meta(pseq_prune),format = 'count')
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="publication",x = "dmn", size = "n", fill = "n")+scale_fill_viridis_b(option = "B")

library(ggforce)
# Do the plotting
a<-as_tibble(my_tbl)%>% group_by(case, dmn) %>%
  summarise(avg = mean(n))
a$dmn<-as_vector(a$dmn)
a$case<-as_vector(a$case)
sam<-meta(pseq_prune)
sam$accesion<-rownames(sam)

sam2<-as_tibble(sam)%>%select(accesion,dmn,case,publication,sample_type )
vec<-as_vector(c("case","sample_type","publication" ))

ggplot(sam2,mapping = aes(x =vec,y = dmn,id=accesion))+
  geom_parallel_sets(aes(fill = dmn), alpha = 0.3, axis.width = 0.1) 
#  geom_parallel_sets_axes(axis.width = 0.1) +
#  geom_parallel_sets_labels(colour = 'white')
###########################################
###Dont forget to save you shit HERE#######
###########################################
save.image(file = "go_terms_dmm.rdata")



count<-abundances(pseq_prune)
#count
#rowMeans(count,na.rm = T)

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