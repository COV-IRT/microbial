library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(tsnemicrobiota)
counts<-read.csv("COUNTS.csv", header = T, row.names = 1, sep = ",")
taxonomy<-read.csv("TAXONOMY.CSV", header = T, row.names = 1, sep = ",")
metadata<-read.csv("metadata.csv", header = T, row.names = 1, sep = ",")


#Make the phyloseq table
counts<-otu_table(counts, taxa_are_rows = T)
taxonomy<-tax_table(as.matrix(taxonomy))
metadata<-sample_data(metadata)
env<-phyloseq(counts, taxonomy, metadata)
#subset the phylo object by trimmomatic type and relabund

env25<-subset_samples(physeq = env, Trim=="trim25")
env25num<-subset_samples(physeq = env25, TrimType=="num")
env25rel<-subset_samples(physeq = env25, TrimType=="frac")
env25relbac<-subset_taxa(physeq = env25rel, Domain=="Bacteria")
env25relbactrim<-core(env25relbac, detection = 0.1/100, prevalence = 50/100, include.lowest = FALSE)


#Plot the Phylum level Bar Plot
plot_taxa_composition(x = env25relbactrim, 
                      sample.sort = "SampleID", 
                      taxonomic.level = "Phylum", 
                      transform = 'compositional', 
                      otu.sort = 'abundance', 
                      plot.type = 'barplot', 
                      average_by = NULL, 
                      x.label = "SampleID",
                      verbose = TRUE) +
  ggtitle(label = 'Environmental Microbiome (Phylum level)', 
          subtitle = "detection = 0.1/100, prevalence = 50/100")


#Plot the top 25 taxa at a species level
env25relbactop<-prune_taxa(x = env25rel,taxa = top_taxa(x = env25relbactrim, n = 25))
plot_taxa_composition(x = env25relbactop, 
                      sample.sort = "SampleID", 
                      taxonomic.level = "Species", 
                      transform = 'compositional', 
                      otu.sort = 'abundance', 
                      plot.type = 'barplot', 
                      average_by = NULL, 
                      x.label = "SampleID",
                      verbose = TRUE) +
  ggtitle(label = 'Environmental Microbiome (Species level)', 
          subtitle = "Top 25 Species")


#t-SNE phyloseq
env25TSNE<-tsne_phyloseq(physeq = env25relbac, 
                         distance = "bray", 
                         perplexity = 50, 
                         dimensions = 2, 
                         precomputed_distance = NULL, 
                         verbose = T, 
                         rng_seed = 314)
sample_data(env25relbac)
p<-plot_tsne_phyloseq(physeq = env25relbac, 
                      tsne_obj = env25TSNE, 
                      color = "SampleID", 
                      title = 't-SNE Bray', justDF = T)

class(p)
title<-"t-SNE (Bray)"
ggplot(p, mapping = aes(X1,X2,color=SampleID, size="0.5"))+
  ggtitle('t-SNE Bray')+
  geom_point()+
  theme_minimal()
                                                                      