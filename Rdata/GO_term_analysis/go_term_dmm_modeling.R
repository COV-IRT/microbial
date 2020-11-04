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
list.files()
#OSF derived combined output from seqscreen
raw<-as_tibble(read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T))

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


####There are multiple processes and values for a single sample so you cant convert the sample to columns
#make individual tibbles for biological processes and molecular fxn
bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")

#make individual tibbles for each type (bac, euk, term, arc, vir, etc)
bio_bac<-bio%>%filter(type=="bac")%>%select(-type)
bio_term<-bio%>%filter(type=="term")%>%select(-type)
mol_bac<-mol%>%filter(type=="bac")%>%select(-type)
mol_term<-mol%>%filter(type=="term")%>%select(-type)

#subselect tibbles for only the counts and go terminology
bio_bac_counts<-bio_bac%>%select(-c(namespace,depth,name))
bio_bac_tax<-bio_bac%>%select(GO_term,namespace,depth,name)

#convert them to dataframes for downstream import to phylsoeq
bio_bac_counts<-data.frame(bio_bac_counts, row.names=1)
bio_bac_tax<-data.frame(bio_bac_tax, row.names=1)


################################################################
# PART 2: MODELING, ANALYSIS, AND VISUALIZATION
################################################################
library(phyloseq)
library(microbiome)
#library(microbiomeutilities)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(ggsci)

#make a phyloseq object
bio_bac_counts<-otu_table(bio_bac_counts, taxa_are_rows = T)
bio_bac_tax<-tax_table(as.matrix(bio_bac_tax), errorIfNULL = T)
bio_bac_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata.txt",header = T, sep = "\t",row.names = 1))
#a little regex to fix the stupid filename
rownames(bio_bac_sam)<-rownames(bio_bac_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")
bio_bac_sam<-sample_data(bio_bac_sam)
bio_bac_physeq<-phyloseq(bio_bac_counts,bio_bac_tax,bio_bac_sam)

#########################################
# Pre-processing
#########################################
#make a copy of the original phyloseq object because you are probably gonna mess it up
pseq<-bio_bac_physeq
# filter out unwanted sampled
pseq<-subset_samples(pseq, sample_type!="neg_control")
pseq<-subset_samples(pseq, sample_type!="Unknown")

#maybe I want to subset for the depth....not sure yet
#pseq_prune<-subset_taxa(physeq = pseq_prune,depth >10)
pseq_prune <- prune_taxa(taxa_sums(pseq) > 10, pseq)
pseq_prune <- prune_samples(sample_sums(pseq_prune) > 10, pseq_prune)
pseq_prune<-tax_glom(physeq = pseq_prune, taxrank = "name")

####################Do NOT EXECUTE THIS CODE#################################################
#write.table(x = rownames(data.frame(otu_table(pseq_prune))), file = "test.tsv",sep = "\t")
#this was causing major issues with the rownames
#taxa_names(pseq_prune)<-get_taxa_unique(pseq_prune,taxonomic.rank = "name" )
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
#convert counts to a matrix
dat <- abundances(pseq_prune)
count <- as.matrix(t(dat))
length(taxa_names(pseq_prune))
#Fit the model
fit <- mclapply(1:10, dmn, count = count, verbose=TRUE)
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
#best <-fit[[4]]
save.image(file = "go_terms_dmm.rdata")

#make a heatmap visualization of the cluster
heatmapdmn(count, fit[[1]], best,ntaxa = 50,
           transform =log2, lblwidth = 0.2 * nrow(count))

#print out the theta values
mixturewt(best)
#save datasheet that show which GO terms contributed to each dmm group
write.table(fitted(best),"GO_TERMS_DMM_contributions.tsv", sep="\t")
#save a datasheet that identifies which sample belongs to which dmm group
ass <- apply(mixture(best), 1, which.max)
write.table(ass,"GO_TERMS_DMM_groups.tsv",sep="")

#make a copy of the phyloseq object so you dont jack it up
physeq<-pseq_prune
meta<-data.frame(sample_data(physeq))
meta$dmm<-ass

# #Go through each dmm cluster and display Go_term contributions to each cluster
# for (k in seq(ncol(fitted(best)))) 
# {
#   d <- melt(fitted(best))
#   colnames(d) <- c("GO", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange GOs by assignment strength
#     arrange(value) %>%
#     mutate(GO = factor(GO, levels = unique(GO))) %>%
#     # Only show the most important drivers
#     filter(abs(value) > quantile(abs(value), 0.8))     
#   
#   p <- ggplot(d, aes(x = GO, y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     labs(title = paste("Top drivers in  : GO  Terms cluster type ", k, sep = ""))
#   #paste(p,k, sep = "")<-p
#   #print(k)
#   print(p)
#   #print(paste(p,k, sep=""))
# }
# 

#add the dmm group to the metadata
sample_data(pseq_prune)$dmn<-ass
#melt the phyloseq object into tidy form
tmp<-psmelt(pseq_prune)
tmp<-as_tibble(tmp)

#subset the dataset to only include the case, Go_term, count, and dmm group.
#obtain the avergage count for each Go term
#order the go terms from hight to lowest count

d2<-tmp %>%
  select(case,name,Abundance, dmn)%>%
  group_by(name,case, dmn) %>%
  summarise(avg = mean(Abundance)) %>%
  arrange(desc(avg))

save.image(file = "go_terms_dmm.rdata")

#move each dmm group into a colum of its own
d3<-tidyr::spread(d2,dmn, avg)

#get the total count of the go terms and oder from greates to lowest
d3[3:6]
d3$tot<-rowSums(d3[3:6], na.rm = T)
d3<-d3%>%arrange(desc(tot))
d3$tot<-NULL

d3<-d3%>%gather(data = d3,avg,3:6)
colnames(d3)<-c("name","case", "dmn","avg")
d4<-d3%>%top_n(100, c(avg))
unique(d4$name)
ggballoonplot(d4, y ="name",x = "case", size = "avg", facet.by = "dmn",fill = "avg")+scale_fill_viridis_c()

