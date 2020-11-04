library(phyloseq)
library(microbiome)
library(ggplot2)
library(microbiomeutilities)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
getwd()
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("sample_metadata_merged.txt", header=T, row.names=1,
                              check.names=F, sep="\t")

# first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(sample_info_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

########################################
# ok now lets filter out the samples and ASVs we dont need
############################################
#remove all the mock community and negative control samples
ASV_physeq<-subset_samples(ASV_physeq, Subject=="Human")
#remove all of the fistula samples that are from Day 0
Aagaard<-subset_samples(ASV_physeq,Dataset=="Aagaard")
Doyle<-subset_samples(ASV_physeq, Dataset=="Doyle")
Aagaard<-subset_samples(physeq = Aagaard,Day_num>0)
Aagaard<-subset_samples(physeq = Aagaard,Day_num!=NA)
ASV_physeq<-merge_phyloseq(Aagaard,Doyle)

#ok lets remove all non bacterial ASVs
ASV_physeq<-subset_taxa(ASV_physeq, Kingdom=="Bacteria")
ASV_physeq
#Now lets aglomerate to a genus level
ASV_physeq_genus<-tax_glom(ASV_physeq, taxrank = "Genus")
ASV_physeq_genus


# ASV_physeq_genus
# ASV_physeq_genus_rel <- microbiome::transform(ASV_physeq_genus, "clr")
# ASV_physeq_genus_rel_core<-core(x = ASV_physeq_genus_rel, detection = 5/100, prevalence = 5/100)
# ASV_physeq_genus_rel_core<-prune_taxa(taxa_sums(ASV_physeq_genus_rel_core) > 0,ASV_physeq_genus_rel_core)
# ASV_physeq_genus_rel_core<-prune_samples(sample_sums(ASV_physeq_genus_rel_core)>0, ASV_physeq_genus_rel_core)
# ASV_physeq_genus_rel_core
# taxa_names(ASV_physeq_genus_rel_core)<-get_taxa_unique(ASV_physeq_genus_rel_core, taxonomic.rank = "Genus", errorIfNULL = T) 


physeq<-ASV_physeq
meta<-sample_data(physeq)
dat <- abundances(physeq)
count <- as.matrix(t(dat))
fit <- lapply(1:3, dmn, count = count, verbose=TRUE)


fit <- mclapply(X = 1:10,FUN =  dmn, count = count,verbose=TRUE)
fit
#Check the model fit with different number of mixture componenets using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)
best <- fit[[which.min(lplc)]]
best
mixturewt(best)
ass <- apply(mixture(best), 1, which.max)
write.table(ass,"Doyle_Aaagaard_GENUS_DMM.tsv",sep="")
write.table(fitted(best),"Doyle_Aagaard_fitted_best_GENUS", sep="")

meta<-data.frame(sample_data(physeq))
meta$dmm<-ass
for (k in seq(ncol(fitted(best)))) 
{
  d <- melt(fitted(best))
  colnames(d) <- c("ASV", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange ASVs by assignment strength
    arrange(value) %>%
    mutate(ASV = factor(ASV, levels = unique(ASV))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = ASV, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers in ",type," : community type ", k, sep = ""))
  #paste(p,k, sep = "")<-p
  #print(k)
  print(p)
  #print(paste(p,k, sep=""))
}

