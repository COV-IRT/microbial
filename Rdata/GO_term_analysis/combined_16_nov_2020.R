library(phyloseq)
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(ggsci)
######################################################
#    PART 1: DATA WRANGLING
######################################################
"###################MIKE LEE PLEASE HELPPPP#############################
########WHY ARE THE GO_TERM NAMES BREAKING THE DATAFRAMES!?!?!#########"

#load libraries
library(phyloseq)
library(tidyverse)

#OSF derived combined output from seqscreen
#raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T)
#dim(raw_df)
# [1] 26966  2020

# this was supposed to have 47,234 lines (or 47,233 rows if read in with header):
# wc -l Combined_BALF_GO_Terms.tsv
#  47234 Combined_BALF_GO_Terms.tsv
# this is due to some having quotes in them that R is trying to interpret as field-quoting characters, e.g.:
# problem will occur in the "name" column, where there can be quotes mixed in, e.g.:
# GO:0006458	biological_process	3	'de novo' protein folding
# or whatever the F this is
# GO:1900549	biological_process	4	N',N'',N'''-triacetylfusarinine C metabolic process
# even though both of those close properly, this one wouldn't
# GO:0061146	biological_process	4	Peyer's patch morphogenesis
# i found those peeking at the command line, e.g.: grep "'" Combined_BALF_GO_Terms.tsv | head
# that is likely causing one or more of these "name" entries to be thousands of rows long

# we can turn that off by setting the quoting characters to nothing, quote=""
#raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T, quote="")
# now we get an error that says
# Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#   line 35939 did not have 2020 elements

# telling us something is up with that row (or the one before or after it depending on how things are being counted)
# peaking there at the command line:
# sed -n '35938,35940p' Combined_BALF_GO_Terms.tsv | cut -f 1-4
# GO:0000982	molecular_function	3	DNA-binding transcription factor activity, RNA polymerase II-specific
# GO:0001133	molecular_function	3	DNA-binding transcription factor activity, RNA polymerase II-specific
# GO:1904067	molecular_function	3	ascr#2 binding

# ah, there is a '#' in that one, nice GO...
# so it's reading that as a comment and then not able to finish that row, we can turn off commenting characters like so:
# raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = "")
# now the whole thing reads in with the expected amount of rows:
# dim(raw_df)

# it's super-valuable to check what we read in matches the expected number of rows and columns :)

# doing it same as you were from here on, just with those added arguments to the initial read.table call


raw<-as_tibble(read.table("Combined_BALF_GO_Terms_parent_propagated.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
raw
# A tibble: 47,233 x 2,020     # good so far now

#do a little regex and fix some stuff
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

#transform the raw table by type of count (euk, term, bac, arc)
df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","Bacteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)

df
####There are multiple processes and values for a single sample so you cant convert the sample to columns
#make individual tibbles for biological processes and molecular fxn
bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")


bio
#make individual tibbles for each type (bac, euk, term, arc, vir, etc)
bio_bac<-bio%>%filter(type=="bac")%>%select(-type)
bio_term<-bio%>%filter(type=="term")%>%select(-type)
mol_bac<-mol%>%filter(type=="bac")%>%select(-type)
mol_term<-mol%>%filter(type=="term")%>%select(-type)

bio_bac
#subselect tibbles for only the counts and go terminology
bio_bac_counts<-bio_bac%>%select(-c(namespace,depth,name))
bio_bac_tax<-bio_bac%>%select(GO_term,namespace,depth,name)
mol_bac_counts<-mol_bac%>%select(-c(namespace,depth,name))
mol_bac_tax<-mol_bac%>%select(GO_term,namespace,depth,name)

#convert them to dataframes for downstream import to phylsoeq
bio_bac_counts<-data.frame(bio_bac_counts, row.names=1)
bio_bac_tax<-data.frame(bio_bac_tax, row.names=1)
mol_bac_counts<-data.frame(mol_bac_counts, row.names=1)
mol_bac_tax<-data.frame(mol_bac_tax, row.names=1)
head(mol_bac_tax)

bio_bac_counts_phy <- otu_table(bio_bac_counts, taxa_are_rows=TRUE)
bio_bac_tax_phy <- tax_table(as.matrix(bio_bac_tax), errorIfNULL=TRUE)
mol_bac_counts_phy<-otu_table(mol_bac_counts, taxa_are_rows = T)
mol_bac_tax_phy<-tax_table(as.matrix(mol_bac_tax), errorIfNULL = T)

# checking nothing similar is happening with the sample info table (though probably not because it has quoting characters)
# number rows (minus one if reading in as with a header)
# wc -l Combined_BALF_GO_Terms_metadata.txt
#      168 Combined_BALF_GO_Terms_metadata.txt
# number of columns
# head -n 1 Combined_BALF_GO_Terms_metadata.txt | tr "\t" "\n" | wc -l
#        70

bio_bac_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata.txt",header = T, sep = "\t",row.names = 1))
# [1] 167  70    # good

#a little regex to fix the stupid filename
rownames(bio_bac_sam)<-rownames(bio_bac_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

# making physeq object

bio_bac_pseq <- phyloseq(bio_bac_counts_phy, bio_bac_tax_phy, sample_data(bio_bac_sam))
mol_bac_pseq<-phyloseq(mol_bac_counts_phy,mol_bac_tax_phy, sample_data(bio_bac_sam))

bac_pseq<-merge_phyloseq(bio_bac_pseq,mol_bac_pseq)

bac_pseq_no_neg<-subset_samples(bac_pseq, sample_type!="neg_control")
bac_pseq_no_neg<-subset_samples(bac_pseq_no_neg, sample_type!="Unknown")


#This seems to be breaking the go term names

bac_pseq_prune <- prune_taxa(taxa_sums(bac_pseq_no_neg) > 10, bac_pseq_no_neg)

#### WARNING, IF YOU DO THIS YOU WILL FREEZE THE CONSOLE AND MAKE IT LAGGY####

bac_pseq_prune <- prune_samples(sample_sums(bac_pseq_prune) > 10, bac_pseq_prune)
#bac_pseq_prune_name<-tax_glom(bac_pseq_prune, taxrank = "name")
bac_pseq_prune

#NAILLLLED IT
bac_pseq_prune_deeper<- subset_taxa(as.numeric(data.frame(tax_table(bac_pseq_prune))$depth)>=10, physeq = bac_pseq_prune) 
bac_pseq_prune_deeper<- prune_taxa(taxa_sums(bac_pseq_prune_deeper) > 10, bac_pseq_prune_deeper)
bac_pseq_prune_deeper<- prune_samples(sample_sums(bac_pseq_prune_deeper) > 10, bac_pseq_prune_deeper)
bac_pseq_prune_deeper # otu_table()   OTU Table:         [ 50 taxa and 71 samples ] not bad

#I think they are using the depth as a character classs rather thana numeric whic is why only 1 and 10 are working for subset
library(tidyverse)

#THIS IS A HACK AND i NEED TO FIND A BETTER WAY TO VECTORIZE THE NAMES INSTEAD OF AGGLOMERATING THE GO TERMS BY NAMES AND LOSING THE RESOLUTION
#OF HAVING BEING ABLE TO SWAP THE GO TERMS WITH THE NAMES IN THE VISUALIZATION WITHOUT LOSING THE RESOLUTION OF AGGLOMERATING GO TERMS WITH THE SAME NAME
#AGLOMERATE BY NAME BECAUSE THERRE ARE ONLY 4 GO TAGS THAT HAVE THE SAME NAME
bac_pseq_prune_deeper2<-tax_glom(bac_pseq_prune_deeper2, taxrank = "name")
#CHANGE THE GO TERM NAME FOR THE IDENTIFIER
taxa_names(bac_pseq_prune_deeper2)<-get_taxa_unique(bac_pseq_prune_deeper2,taxonomic.rank = "name" )

############################################################################################

#########################################
# Dataset summary visulizations
#########################################

# mypal3<-c("#3B4992FF", "#BB0021FF", "#5F559BFF", "#CC0C00FF", "#EE0000FF","#008B45FF",  "#631879FF", "#008280FF", "#5C88DAFF", "#84BD00FF", "#FFCD00FF", "#7C878EFF", "#00B5E2FF", "#00AF66FF")
# p <- plot_frequencies(sample_data(pseq), "publication","sample_type")+scale_fill_manual(values = mypal3)
# q <- plot_frequencies(sample_data(pseq), "publication","case")+scale_fill_manual(values =c( "#008B45FF","#3B4992FF" ,"#EE0000FF"))#green-blue-red 

#########################################
# DMM modeling
#########################################
library(microbiome)
library(DirichletMultinomial)

#convert counts to a matrix
dat <- abundances(bac_pseq_prune_deeper2)
count <- as.matrix(t(dat))
length(taxa_names(bac_pseq_prune_deeper2))
bac_pseq_prune_deeper2 #854 TAXA BY 133 SAMPLES TOTALLY DOABLE ON A LOCAL IN A REASONABLE AMT OF TIME
#Fit the model
fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)
fit

#Check the model fit with different number of mixture componenets using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

#identify the number of clusters that best fits the model
best <- fit[[which.min(lplc)]]
best
best <-fit[[4]]
save.image(file = "bac_go_terms_dmm.rdata")

#make a heatmap visualization of the cluster
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =log2, lblwidth = 0.2 * nrow(count))
#make a heatmap visualization of the cluster
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =sqrt, lblwidth = 0.2 * nrow(count))
#print out the theta values
mixturewt(best)
#save datasheet that show which GO terms contributed to each dmm group
write.table(fitted(best),"combined_bac_GO_TERMS_DMM_contributions.tsv", sep="\t")
#save a datasheet that identifies which sample belongs to which dmm group
ass <- apply(mixture(best), 1, which.max)
write.table(ass,"combined_bac_GO_TERMS_DMM_groups.tsv",sep="")

#make a copy of the phyloseq object so you dont jack it up
physeq<-bac_pseq_prune_deeper2
meta<-data.frame(sample_data(physeq))
meta$dmm<-ass
library(reshape2)
#Go through each dmm cluster and display Go_term contributions to each cluster
for (k in seq(ncol(fitted(best))))
{
  d <- melt(fitted(best))
  colnames(d) <- c("GO", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange GOs by assignment strength
    arrange(value) %>%
    mutate(GO = factor(GO, levels = unique(GO))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.99))
  
  p <- ggplot(d, aes(x = GO, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers in  : GO  Terms cluster type ", k, sep = ""))+
           ggsave(filename = paste("contibution_bar_plot_",k,sep = ""),device = "png")
  #paste(p,k, sep = "")<-p
  #print(k)
  print(p)
  #print(paste(p,k, sep=""))
}


#add the dmm group to the metadata
# sample_data(bac_pseq_prune_deeper2)$dmn<-ass
# #melt the phyloseq object into tidy form
#  tmp<-psmelt(bac_pseq_prune_deeper2)
#  tmp<-as_tibble(tmp)
# 
# #subset the dataset to only include the case, Go_term, count, and dmm group.
# #obtain the avergage count for each Go term
# #order the go terms from hight to lowest count
# 
# d2<-tmp %>%
#   select(case,name,Abundance, dmn)%>%
#   group_by(name,case, dmn) %>%
#   summarise(avg = mean(Abundance)) %>%
#   arrange(desc(avg))
# 
save.image(file = "combined_bac_go_terms_dmm.rdata")
# 
# #move each dmm group into a colum of its own
# d3<-tidyr::spread(d2,dmn, avg)
# 
# #get the total count of the go terms and oder from greates to lowest
# d3[3:5]
# d3$tot<-rowSums(d3[3:5], na.rm = T)
# d3<-d3%>%arrange(desc(tot))
# d3$tot<-NULL
# 
# d3<-d3%>%gather(data = d3,avg,3:5)
# colnames(d3)<-c("name","case", "dmn","avg")
# d4<-d3%>%top_n(100, c(avg))
# unique(d4$name)
# library(ggpubr)
# ggballoonplot(d4, y ="name",x = "case", size = "avg", facet.by = "dmn",fill = "avg")+scale_fill_viridis_c()
# 
#bac_pseq_prune_deep
tmp<-psmelt(bac_pseq_prune_deeper2)
tmp<-as_tibble(tmp)

tmp
d2<-tmp %>%
  select(case,name,Abundance, dmn)%>%
  group_by(name,case, dmn) %>%
  summarise(avg = mean(Abundance),std = sd(Abundance)) %>%
  group_by(name, case,dmn, std)%>%
  arrange(desc(std))
d2
#d3<-tidyr::spread(d2,dmn, avg)
d3<-tidyr::spread(d2,case,avg)

d3[4:6]
#library(mosaic)
#d3$tot<-sd(d3[3:5])
#d3$tot<-rowSums(d3[3:5], na.rm = T)
d3<-d3%>%arrange(desc(std))
rank<-as.character(d3$name)
dim(d3)
d4<-d3[1:800,]

d4<-d4%>%
  gather(data = d4,avg,4:6)

colnames(d4)<-c("name","dmn","std", "group","avg")

#d4<-d3%>%group_by(name,group)
d4<-d4%>%filter(!is.na(avg))
library(ggpubr)
d4
ggballoonplot(data = d4, y ="name",x = "dmn", size = "avg", facet.by = "group",fill = "avg")+scale_fill_viridis_c()
write.table(meta(bac_pseq_prune_deeper2),file = "combined_pseq_meta.tsv",sep = '\t')
library(mosaic)

meta(pseq_prune)
head(tmp$OTU)
head(tmp$Abundance)
head(tmp$dmn)



my_tbl<-tally(~dmn,data =meta(bac_pseq_prune_deeper2),format = 'count')
my_tbl
library(microbiome)
my_tbl<-tally(~case+ publication+dmn,data =meta(bac_pseq_prune_deeper2),format = 'count')
my_tbl
a<-as_tibble(my_tbl)

ggballoonplot(data = a, y ="case",x = "dmn", facet.by = "publication",size = "n", fill = "n")+
  scale_fill_viridis_c(option = "C")+
  ggsave(filename = "pub_vs_disease_vs_dmm_dotplot.png",
         device = "png",
         #width = "8",
         #height = "6",
         #units = "in",
         dpi = 600,
         path = "D:/github/microbial/Rdata/GO_term_analysis/Figures/")

my_tbl<-tally(~case+ dmn,data =meta(bac_pseq_prune_deeper2),format = 'count')
my_tbl
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="case",x = "dmn", size = "n", fill = "n")+
  scale_fill_viridis_c(option = "B")+
  ggsave(filename = "disease_vs_dmm_dotplot.png",
         device = "png",
         #width = "8",
         #height = "6",
         #units = "in",
         dpi = 600,
         path = "D:/github/microbial/Rdata/GO_term_analysis/Figures/")

my_tbl<-tally(publication ~ dmn,data =meta(bac_pseq_prune_deep),format = 'count')
my_tbl
a<-as_tibble(my_tbl)
ggballoonplot(data = a, y ="publication",x = "dmn", size = "n", fill = "n")+scale_fill_viridis_b(option = "B")

###########################################
###Dont forget to save you shit HERE#######
###########################################
save.image(file = "bac_go_terms_dmm.rdata")


library(dplyr)
count<-abundances(bac_pseq_prune_deeper2)
select <- order(rowMeans(count),decreasing=TRUE)
select2<-sqrt((count)[select,])
tmp<-rownames(select2)
#select2$mol<-tmp
select3<-as.table(select2)
test<-t.test(select3)
test
library(matrixTests)
library(genefilter)
select3<-as_tibble(select2)%>%summarise(std=rowFtests(select2))%>% arrange(desc(std))
select3
select2
select2$mol<-tmp


select2
select2<-as_tibble(select2)
select2
select2$mean<-rowMeans(select2)
select3
select2<-filter(.data = select2,rowMeans(select2)<10)

sam<-data.frame(sample_data(bac_pseq_prune_deeper2))
df<-as.data.frame(sample_data(bac_pseq_prune_deeper2))
df<-as_tibble(df)
df<-df%>%select(dmn, publication, sample_type,case)#,dmn,body_site)
df
df<-as.data.frame(df)
df

select2<-sqrt((count)[select,])
dim(select2)
colnames(select2) <- colnames(otu_table(bac_pseq_prune_deep))
length(rownames(select2))
length(row.names(df))
row.names(df) <- colnames(select2)


library(ggsci)
mypal <- pal_aaas("default", alpha = 1)(10)
mypal
library("scales")
library(RColorBrewer)
library(viridis)
library(pheatmap)
df_row<-as.data.frame(fitted(best))
df<-data.frame(df)
df
colnames(df)<-c("dmm_cluster", "Publication","Sample_Type","Case")
# Specify colors
ann_colors = list(
  dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF","5"="#008280FF"),
  Publication=c("Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
  Sample_Type =c("COVID_19"="firebrick","Healthy"="forestgreen","Community_acquired_pneumonia"="dodgerblue4",
                 "Obese_Asthma"="goldenrod3",
                 "Obese_Smoker"="goldenrod4", 
                 "Obese"="goldenrod1",
                 "Asthma"= "darkorange2", 
                 "Asthma_Smoker"="darkorange4",
                 "Asthma_Ex_smoker"="darkorange3",
                 "Smoker"="gray27",
                 "Obese_Asthma_Smoker"="black"),
  Case=c("COVID19"="firebrick", "Control_Healthy"="forestgreen","Control_Sick"="dodgerblue4"))


xx <- pheatmap(mat = select2,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),
               annotation_col=df,
               annotation_colors = ann_colors,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               annotation_row = df_row)
xx
# pheatmap(mat, 
#          color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
#          kmeans_k = 5, 
#          breaks = NA, 
#          border_color = "grey60",
#          cellwidth = NA, 
#          cellheight = NA, 
#          scale = "none", 
#          cluster_rows = TRUE,
#          cluster_cols = TRUE, 
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean", 
#          clustering_method = "complete",
#          clustering_callback = identity2, 
#          cutree_rows = NA, 
#          cutree_cols = NA,
#          treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
#                                  50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
#                                                                    cluster_cols, 50, 0), 
#          legend = TRUE,
#          legend_breaks = NA,
#          legend_labels = NA,
#          annotation_row = NA, 
#          annotation_col = NA,
#          annotation = NA, 
#          annotation_colors = NA, 
#          annotation_legend = TRUE,
#          annotation_names_row = TRUE, 
#          annotation_names_col = TRUE,
#          drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
#          fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
#          angle_col = c("270", "0", "45", "90", "315"), display_numbers = F,
#          number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8
#          * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
#          labels_col = NULL, filename = NA, width = NA, height = NA,
#          silent = FALSE, na_col = "#DDDDDD", ...)
