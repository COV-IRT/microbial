library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

#load the phyloseq object that has the decontamed samples and no neg controls
pseq_decontam_no_neg_core
deseq_counts <- phyloseq_to_deseq2(physeq = pseq_decontam_no_neg_core,design = ~1)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
#Do a VST transformation on the deseq_counts
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)

#make a distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$case[order.dendrogram(euc_dend)])
labels(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_tax_phy <- tax_table(pseq_decontam_no_neg_core)
sample_info_tab_phy <- sample_data(pseq_decontam_no_neg_core)

vst_physeq <- phyloseq(vst_count_phy,vst_tax_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="case") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=publication, hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")  
#  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
count_tab<-data.frame(otu_table(pseq_decontam_no_neg_core))
rarecurve(t(count_tab), step=100, lwd=2, ylab="Read Counts", label=F)

# and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))


#################################
#Richness and diversity estimates
#################################

# first we need to create a phyloseq object using our un-transformed count table


# and now we can call the plot_richness() function on our phyloseq object
plot_richness(pseq_decontam_no_neg_core, color="case", measures=c("Chao1", "Shannon")) + 
  theme(legend.title = element_blank())

plot_richness(pseq_decontam_no_neg_core, x="case", color="case", measures=c("Chao1", "Shannon")) + 
  theme(legend.title = element_blank())

sam<-sample_data()
anova(betadisper(euc_dist, sample_info_tab$case)) # 0.002
betadisper(euc_dist, sample_info_tab$case)
