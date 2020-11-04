library(phyloseq)
setwd("d:/github/microbial/metaphlan3_launcher/profiles/")

counts<-read.table("counts.txt",header = T, row.names = 1, sep = "\t")
tax<-read.table("tax.txt",header = T, row.names = 1, sep = "\t")
sam<-read.table("metaphlan3_sample_metadata.txt",header = T, row.names = 1, sep = "\t")
tax<-as.matrix(tax)
counts<-otu_table(counts,taxa_are_rows = T)
tax<-tax_table(tax,errorIfNULL = T)
sam<-sample_data(sam)
ps<-phyloseq(counts,tax,sam)
ps<-subset_taxa(physeq = ps,Kingdom=="k__Bacteria")
ps<-subset_samples(physeq = ps, case!="NA")
ps<-subset_samples(physeq = ps, case!="Control_Unknown")
ps<-subset_samples(physeq = ps, tissue=="BALF")
ps_phylum<-tax_glom(ps,taxrank = "Phylum")
taxa_names(ps_phylum)<-get_taxa_unique(physeq = ps_phylum,taxonomic.rank = "Phylum")
plot_composition(ps_phylum, otu.sort = "abundance",group_by =  "case",verbose = T)
ps_genus<-tax_glom(physeq = ps, taxrank = "Species",NArm = T)
ps_genus_comp<-prune_samples(x = ps_genus,samples = sample_sums(ps_genus)>0)

taxa_names(ps_genus)<-get_taxa_unique(physeq = ps_genus,taxonomic.rank = "Species")
sample_sums(ps_genus)
ps_genus
ps_genus_comp<-prune_samples(x = ps_genus,samples = sample_sums(ps_genus)>0)
ps_genus_comp
ps_genus_comp<-subset_taxa(physeq = ps_genus_comp,taxa_sums(ps_genus_comp)>0)
ps_genus_comp
ps_core<-core(x = ps_genus_comp, detection = 5/100,prevalence = 5/100)
ps_core
library(microbiome)
library(ggsci)

my_pal<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000',
          '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

plot_composition(x = ps_core, otu.sort = "abundance",group_by = "case")+scale_fill_manual(values = my_pal)

