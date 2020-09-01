library("phyloseq")
library("microbiome")
library("DESeq2")
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

#I decided to remove these samples from the analysis because they had a strong batch effect
#and were outliers
ps_balf<-subset_samples(physeq = ps_balf, publication!="Huang et al. 2019")
#lets preprocess a bit and find out if there are any taxa that arent present in any samples
any(taxa_sums(ps_balf) == 0) #TRUE
ps_balf <- prune_taxa(taxa_sums(ps_balf) > 0, ps_balf)#lets get rid of those
#How many taxa are there?
ntaxa(ps_balf)#5328

######################################################
#Ok, lets take a look at the sequencing depth###
################################################
SeqDepth = colSums(otu_table(ps_balf))
sample_data(ps_balf)$SeqDepth = SeqDepth
qplot(log10(SeqDepth), geom = "histogram") + theme_bw()
library(ggpubr)
sam<-data.frame(sample_data(ps_balf))
sam$publication<-as.factor(sam$publication)
sam$SeqDepthlog10<-log10(SeqDepth)
set.seed(1234)
wdata = data.frame(
  sex = factor(rep(c("F", "M"), each=200)),
  weight = c(rnorm(200, 55), rnorm(200, 58)))

gghistogram(data = sam,
            x = "SeqDepthlog10",
            fill ="publication", 
            title = "Sequencing Depth",
            add_density = TRUE)+  scale_fill_brewer(palette = "Set1",direction = -1)+theme_bw()

gghistogram(data = sam,
            x = "SeqDepthlog10",
            fill ="case", 
            title = "Sequencing Depth by case",
            add_density = TRUE)+  scale_fill_brewer(palette = "Set1")+theme_bw()



######VST tranformation############

deseq_counts <- phyloseq_to_deseq2(ps_balf,design = ~case)
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
dend_cols <- as.character(sam_tab$publication[order.dendrogram(euc_dend)])
labels(euc_dend) <- dend_cols
plot(euc_dend, ylab="VST Euc. dist.")

plot(euc_clust)
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sam_tab$case[order.dendrogram(euc_dend)])

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

plot_ordination(vst_physeq, vst_pcoa, color="case", shape="case") + 
  scale_shape_manual(values=1:nlevels(sam_tab$publication)) +
  geom_point(size=2) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")+ 
  theme_bw()+scale_color_brewer(palette = "Set1")
a
a+facet_wrap(facets = ~publication)+stat_ellipse(aes(color=case, group=case))


a<-plot_ordination(vst_physeq, vst_pcoa, color="publication", shape="case") + 
  scale_shape_manual(values=1:nlevels(sam_tab$publication)) +
  geom_point(size=2) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA")+ 
  theme_bw()+scale_color_brewer(palette = "Set1")

###Lets try a tsne plot
sample_data(ps_balf)$case<-factor(x = sample_data(ps_balf)$case,levels = c("Control_Healthy","Control_Sick","COVID19"))
library(tsnemicrobiota)
tsne_balf<-tsne_phyloseq(physeq = ps_balf,distance = 'bray',perplexity = 25,verbose = T,rng_seed = 314)
plot_tsne_phyloseq(physeq = ps_balf,
                   tsne_obj = tsne_balf,
                   color = "case",
                   shape = "publication",
                   title = "Bray curtis T-SNE plot by case")+
  stat_ellipse(aes(color=case, group=case))+
  scale_color_brewer(palette = "Set1")+
  theme_bw()
plot_tsne_phyloseq(physeq = ps_balf,
                   tsne_obj = tsne_balf,
                   color = "publication",
                   shape = "case",
                   title = "Bray curtis T-SNE plot by publication")+
  stat_ellipse(aes(color=publication, group=publication))+
  scale_color_brewer(palette = "Set1",direction = -1)+
  theme_bw()


###################Differentially abundant Taxa############
ps_balf
ps_balf_genus<-tax_glom(physeq = ps_balf,taxrank="genus",NArm = T)

taxa_names(ps_balf_genus)<-get_taxa_unique(physeq =ps_balf_genus, taxonomic.rank = "genus")

deseq<-phyloseq_to_deseq2(ps_balf, design = ~case)
deseq <- DESeq(deseq)
resultsNames(deseq)
# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Healthy_vs_COVID19<- results(deseq, alpha=0.01, contrast=c("case", "Control_Healthy","COVID19"))
deseq_res_Sick_vs_COVID19<- results(deseq, alpha=0.01, contrast=c("case", "Control_Sick", "COVID19"))
deseq_res_Healthy_vs_Sick<- results(deseq, alpha=0.01, contrast=c("case", "Control_Healthy","Control_Sick"))
summary(deseq_res_Healthy_vs_COVID19)
summary(deseq_res_Sick_vs_COVID19)
summary(deseq_res_Healthy_vs_Sick)
#subset it for only statistically significant
sigtab_res_Healthy_vs_COVID19 <- deseq_res_Healthy_vs_COVID19[which(deseq_res_Healthy_vs_COVID19$padj < 0.01), ]
sigtab_res_Sick_vs_COVID19 <- deseq_res_Sick_vs_COVID19[which(deseq_res_Sick_vs_COVID19$padj < 0.01), ]
sigtab_res_Healthy_vs_Sick <- deseq_res_Healthy_vs_Sick[which(deseq_res_Healthy_vs_Sick$padj < 0.01), ]

#subset it for only log fold changes >2
sigtab_res_Healthy_vs_COVID19 <- sigtab_res_Healthy_vs_COVID19[which(abs(sigtab_res_Healthy_vs_COVID19$log2FoldChange) >10), ]
sigtab_res_Sick_vs_COVID19 <- sigtab_res_Sick_vs_COVID19[which(abs(sigtab_res_Sick_vs_COVID19$log2FoldChange) >10), ]
sigtab_res_Healthy_vs_Sick <- sigtab_res_Healthy_vs_Sick[which(abs(sigtab_res_Healthy_vs_Sick$log2FoldChange) >10), ]
summary(sigtab_res_Healthy_vs_COVID19)
summary(sigtab_res_Sick_vs_COVID19)
summary(sigtab_res_Healthy_vs_Sick)


sigtab_res_Healthy_vs_COVID19_with_tax <- cbind(as(sigtab_res_Healthy_vs_COVID19, "data.frame"), as(tax_table(ps_balf)[row.names(sigtab_res_Healthy_vs_COVID19), ], "matrix"))
sigtab_res_Sick_vs_COVID19_with_tax <- cbind(as(sigtab_res_Sick_vs_COVID19, "data.frame"), as(tax_table(ps_balf)[row.names(sigtab_res_Sick_vs_COVID19), ], "matrix"))
sigtab_res_Healthy_vs_Sick_with_tax <- cbind(as(sigtab_res_Healthy_vs_Sick, "data.frame"), as(tax_table(ps_balf)[row.names(sigtab_res_Healthy_vs_Sick), ], "matrix"))


# sort it by > change in base_mean
sigtab_res_Healthy_vs_COVID19_with_tax<-sigtab_res_Healthy_vs_COVID19_with_tax[order(sigtab_res_Healthy_vs_COVID19_with_tax$baseMean, decreasing=T), ]
sigtab_res_Sick_vs_COVID19_with_tax<-sigtab_res_Sick_vs_COVID19_with_tax[order(sigtab_res_Sick_vs_COVID19_with_tax$baseMean, decreasing=T), ]
sigtab_res_Healthy_vs_Sick_with_tax<-sigtab_res_Healthy_vs_Sick_with_tax[order(sigtab_res_Healthy_vs_Sick_with_tax$baseMean, decreasing=T), ]

dim(sigtab_res_Healthy_vs_COVID19_with_tax)
dim(sigtab_res_Sick_vs_COVID19_with_tax)
dim(sigtab_res_Healthy_vs_Sick_with_tax)

##########boxplots comparions of DESEQ2 identifies taxa########
library(dplyr)
#merge the dataframes together (this works but needs to be cleaned up)
DeSEQ_tax<-merge.data.frame(sigtab_res_Healthy_vs_COVID19_with_tax,sigtab_res_Sick_vs_COVID19_with_tax, by=0, all=T)
rownames(DeSEQ_tax)<-DeSEQ_tax$Row.names
DeSEQ_tax<-merge(DeSEQ_tax,sigtab_res_Healthy_vs_Sick_with_tax,by=0, all=T)
rownames(DeSEQ_tax)<-DeSEQ_tax$Row.names
DeSEQ_tax$Row.names<-NULL
list<-as_tibble(DeSEQ_tax)

list<-as_tibble(list)%>%distinct(Row.names, .keep_all = F)%>%filter(!is.na(Row.names))
print(list)

rownames(list)<-list$Row.names
list2<-as.character(list$Row.names)
list
ps_balf_comp<- transform(ps_balf, 'compositional')
select<-prune_taxa(x = ps_balf_comp,taxa=list2)
sample_sums(select)
select<-psmelt(select)

select

my_comparisons<-list(c("Control_Healthy","Control_Sick"),c("Control_Sick","COVID19"),c("Control_Healthy","COVID19"),)

ggviolin(data = select, x = "case", y = "Abundance",order = c("Control_Healthy","Control_Sick","COVID19"),
          color = "case",add = "jitter", shape = NULL,title = "Differentially Abundant Taxa",ylab = "Relative Abundance (log10)")+ 
  stat_compare_means(comparisons = my_comparisons, method = "sdfawilcox.test")+
  facet_wrap(nrow = 3,ncol = 13,facets = ~genus)+
  scale_color_brewer(palette = "Set1", direction = -1)+rotate_x_text()+scale_y_log10()

#############Communitity compoistions#############
plot_core(x = ps_balf_comp,plot.type = 'heatmap', prevalences=seq(0.1, 1, .1), detections=seq(0.01, 1, length = 10))
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown","dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown","dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
#subset the top 25 taxa and change the tax_ids to the names
ps_balf_comp_top<-prune_taxa(taxa = top_taxa(x = ps_balf_comp,n = 25),x = ps_balf_comp)
taxa_names(ps_balf_comp_top)<-get_taxa_unique(ps_balf_comp_top,taxonomic.rank = "species")

plot_composition(x = ps_balf_comp_top,
                 sample.sort = 'neatmap',
                 otu.sort = "abundance",
                 x.label = "publication",
                 plot.type = 'barplot',
                 verbose = T,
                 group_by = "case")+
  scale_fill_manual(values = c25)+
  theme_bw()+rotate_x_text()






##############################BETA DIVERSITY#################################################

within_group_difs_anova<-anova(betadisper(euc_dist, meta(vst_physeq)$case)) # 0.002
within_group_difs_permutest<-permutest(betadisper(euc_dist, meta(vst_physeq)$case))
TukeyHSD(betadisper(euc_dist, meta(vst_physeq)$case))


##########################################################
##################COVID19 vs Healthy###################################
##########################################################
ps_balf_not_sick<-subset_samples(physeq = ps_balf, case!="Control_Sick")

deseq_counts <- phyloseq_to_deseq2(ps_balf_not_sick,design = ~case)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 

# running betadisper on just these based on level of alteration as shown in the images above:
permutest(betadisper(euc_dist,sample_data(ps_balf_not_sick)$case))
write.table(x = anova((betadisper(euc_dist,sample_data(ps_balf_not_sick)$case))),file = "COVID_19_vs_Healthy_betadisper.tsv",sep = "\t",append = T)# 0.7

a<-adonis(euc_dist~sample_data(ps_balf_not_sick)$case,strata = sample_data(ps_balf_not_sick)$publication, permutations = 9999) # 0.003
a
capture.output(x = a,file = "COVID19_vs_Healthy_strata_publication_adonis.tsv",append = F)# 0.7
