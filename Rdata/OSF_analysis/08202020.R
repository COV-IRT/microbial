library(phyloseq)
library(microbiome)
library(DESeq2)
library("vegan")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

counts_tab<-read.csv(file = "num_counts.txt",header = T,row.names = 1,sep = "\t")
tax_tab<-as.matrix(read.csv(file = "taxonomy.txt",header = T,row.names = 1,sep = "\t"))
sam_tab<-read.csv(file = "sample_metadata.txt",header = T,row.names = 1,sep = "\t")
counts<-otu_table(counts_tab,taxa_are_rows = T)
tax<-tax_table(tax_tab,errorIfNULL = T)
sam<-sample_data(sam_tab)
ps<-phyloseq(counts,tax,sam)
rank.names(ps)
ps<-subset_taxa(ps, domain=="Bacteria")
ps_balf<-subset_samples(physeq = ps, tissue=="BALF")
ps_pbmc<-subset_samples(physeq = ps, tissue=="PBMC")



deseq_counts <- phyloseq_to_deseq2(ps_balf,design = ~case)
#deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sam_tab$bioproject[order.dendrogram(euc_dend)])
labels(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")


# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sam_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
sam_tab$publication <- factor(sam_tab$publication)

a<-plot_ordination(vst_physeq, vst_pcoa, color="case", shape="case") + 
  scale_shape_manual(values=1:nlevels(sam_tab$publication)) +
  geom_point(size=2) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")+ 
  theme_bw()+scale_color_brewer(palette = "Set1")

a+facet_wrap(facets = ~publication)+stat_ellipse(aes(color=case, group=case))
