## GO term level differential expression analysis using DESeq2

## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggrepel)

input_dir <- tryCatch(
  {
    paste(dirname(rstudioapi::getSourceEditorContext()$path), "/", sep="")
  }, error = function(e) {
    "./"
  }
)
  
setwd(input_dir)

output_dir <- paste(input_dir, "results/", sep="")
dir.create(file.path(output_dir), showWarnings = FALSE)

# Read data file
# Replace all ' and # before reading! XXX Make it automatic.
data_file <- paste(input_dir, "Transformed_BALF_GO_Terms_Subset.tsv", sep="")
noquote(paste("Reading data from", data_file))
data <- read.table(data_file, sep="\t", header=TRUE, row.names=1)
#data[which(data[,1]=="GO:0008663"),]
#which(rownames(data)=="GO:0008663")
#data["GO:0008663",3]
#data[13736,3]

## Columns 1-3 are metadata for the GO term, move to separate table
meta_go_term <- data[1:3]
data <- data[-c(1:3)]
  
## Read metadata table
meta_file <- paste(input_dir, "Metadata.tsv", sep="")
noquote(paste("Reading data from", meta_file))
meta <- read.table(meta_file, sep="\t", header=T, row.names=1)

## RNA-seq count distribution
#ggplot(data)+
#  geom_histogram(aes(x = CRR125934_term_counts), stat = "bin", bins = 200) +
#  xlab("Raw expression counts") +
#  ylab("Number of GO terms")
  
## Play with xlim to see distribution of GO terms
#ggplot(data)+
#  geom_histogram(aes(x = CRR125934_term_counts), stat = "bin", bins = 200) +
#  xlim(-5, 200)  +
#  xlab("Raw expression counts") +
#  ylab("Number of GO terms")
## ^ Low number of counts associated with a large proportion of GO terms

## Modeling count data
## Play with data range
## -- data for all data
## -- data[, 1:2] or data[, 3:4] for replicates 
## -- data[, 1:16] or data[, 17-41] for groups
#mean_counts <- apply(data[, 17:41], 1, mean)
#variance_counts <- apply(data[, 17:41], 1, var)
#df <- data.frame(mean_counts, variance_counts)
#ggplot(df) +
#  geom_point(aes(x=mean_counts, y=variance_counts)) + 
#  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
#  scale_y_log10() +
#  scale_x_log10()
## ^ In all cases mean < variance, so Negative Binomal distribution is a best fit
  
  
## Unify sample type for all CAP
noquote("Unifying sample_type for all CAP")
meta$sample_type[gsub("_.*", "", meta$sample_type) == "CAP"] <- "CAP"

## Transform colnames in data and rownames in meta to match
noquote("Unifying column names in data and row names in metadata")
colnames(data) <- gsub("_.*","", colnames(data))
new_rownames <- c()
for (rn in rownames(meta)) {
  accession <- gsub("\\..*", "", rn)
  sample_name <- gsub(".*\\.","", rn)
  if ( accession %in% colnames(data)) {
    new_rownames <- c(new_rownames, accession)
  }
  else if (sample_name %in% colnames(data)) {
    new_rownames <- c(new_rownames, sample_name)
  }
  else { 
    noquote(paste(rn, "not found in the metadata table.")) 
    new_rownames <- c(new_rownames, rn)
  }
}
rownames(meta) <- new_rownames
  
data <- data[, match(rownames(meta), colnames(data))]

### Check that sample names match in both files
noquote("Checking if all column names in data have corresponding rows in metadata")
if (! all(colnames(data) %in% rownames(meta))) {
  stop("Sample names not matched in data and metadata files")
}

if (! all(colnames(data) == rownames(meta))) {
  stop("Sample names not matched in data and metadata files")
}

## Create DESeq2Dataset object
#dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~sample_type)
noquote("Creating DeSeq matrix")
dds <- DESeqDataSetFromMatrix(countData=data, colData=meta, design=~sample_type)

## Save normalized counts for later use
dds_esf <- estimateSizeFactors(dds)
#sizeFactors(dds)
normalized_counts <- counts(dds_esf, normalized=TRUE)
normalize_count_file <- paste(output_dir, "normalized_counts.tsv", sep="")
noquote(paste("Saving normalized counts to", normalize_count_file))
write.table(normalized_counts, file=normalize_count_file, sep="\t", quote=F, col.names=NA)
  
## Run analysis
noquote("Running analysis")
dds <- DESeq(dds)

## Exploration: uncomment lines to explore dds features 
## Check the size factors
#sizeFactors(dds)
## Total number of raw counts per sample
#colSums(counts(dds))
## Total number of normalized counts per sample
#colSums(counts(dds, normalized=T))
## Plot dispersion estimates
#plotDispEsts(dds)

## COVID19 vs Healthy
## Define contrasts, extract results table, and shrink the log2 fold changes
contrast_covid19 <- c("sample_type", "COVID19_BALF", "Healthy_BALF")
res_covid19_unshrunken <- results(dds, contrast=contrast_covid19, alpha = 0.05)
res_covid19 <- lfcShrink(dds, contrast=contrast_covid19, res=res_covid19_unshrunken, type="ashr")

## Select top 20 GO terms by padj values, add metadata for GO terms and write results
top20_terms_covid19 <- res_covid19 %>% data.frame() %>% arrange(padj) %>% head(n=20)
top20_terms_covid19 <- merge(meta_go_term, top20_terms_covid19, by="row.names") 
top20_terms_covid19 <- arrange(top20_terms_covid19, padj)
row.names(top20_terms_covid19)<- top20_terms_covid19$Row.names
top20_terms_covid19 <- top20_terms_covid19[-1]
go_term_meta_top20_file = paste(output_dir, "top20_go_terms_covid19_healthy.tsv", sep="")
noquote(paste("Saving top20 significant GO terms in COVID19 vs Healthy to", go_term_meta_top20_file))
write.table(top20_terms_covid19, file=go_term_meta_top20_file, sep="\t", quote=F, col.names=NA)

## COVID19 vs CAP
## Define contrasts, extract results table, and shrink the log2 fold changes
contrast_covid19 <- c("sample_type", "COVID19_BALF", "CAP")
res_covid19_unshrunken <- results(dds, contrast=contrast_covid19, alpha = 0.05)
res_covid19 <- lfcShrink(dds, contrast=contrast_covid19, res=res_covid19_unshrunken, type="ashr")

## Select top 20 GO terms by padj values, add metadata for GO terms and write results
top20_terms_covid19 <- res_covid19 %>% data.frame() %>% arrange(padj) %>% head(n=20)
top20_terms_covid19 <- merge(meta_go_term, top20_terms_covid19, by="row.names")
top20_terms_covid19 <- arrange(top20_terms_covid19, padj)
row.names(top20_terms_covid19)<- top20_terms_covid19$Row.names
top20_terms_covid19 <- top20_terms_covid19[-1]
go_term_meta_top20_file = paste(output_dir, "top20_go_terms_covid19_cap.tsv", sep="")
noquote(paste("Saving top20 significant GO terms in COVID19 vs CAP to", go_term_meta_top20_file))
write.table(top20_terms_covid19, file=go_term_meta_top20_file, sep="\t", quote=F, col.names=NA)


## What are common terms in top20_go_terms_covid_healthy.tsv and top20_go_terms_covid_cap.tsv
covid19_healthy_file <- paste(output_dir, "top20_go_terms_covid19_healthy.tsv", sep="")
covid19_healthy_terms <- read.table(covid19_healthy_file, sep="\t", header=T, row.names=1)
covid19_cap_file <- paste(output_dir, "top20_go_terms_covid19_cap.tsv", sep="")
covid19_cap_terms <- read.table(covid19_cap_file, sep="\t", header=T, row.names=1)
common_terms <- merge(covid19_healthy_terms, covid19_cap_terms, by="row.names")[,1]
common_terms <- meta_go_term[common_terms,]
common_terms_file = paste(output_dir, "top_common_covid19_healthy_and_covid19_cap.tsv", sep="")
noquote(paste("Saving common significant GO terms to", common_terms_file))
write.table(common_terms, file=common_terms_file, sep="\t", quote=F, col.names=NA)

# Plot expression for single gene
noquote("Plotting common significant GO terms")
for (term in rownames(common_terms)) {
  plot_data <- plotCounts(dds, term, intgroup="sample_type", returnData = TRUE)
  plot_file = paste(output_dir, gsub(":", "-", term), ".png", sep="")
  png(file=plot_file, width=1200, height=900, res=120)
  plot_to_save <- ggplot(plot_data, aes(x=sample_type, y=count, color=sample_type)) + 
    ggtitle(term) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("normalized count") +
    geom_jitter(size=6) + 
    scale_y_log10()
  print(plot_to_save)
  dev.off()
}

noquote(paste("Plots for", paste(rownames(common_terms), collapse=", "), "saved in", output_dir))


