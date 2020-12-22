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
data_file <- paste(input_dir, "Transformed_BALF_GO_Terms_Subset.tsv", sep="")
noquote(paste("Reading data from", data_file))
file_rows = strtoi(sub("\\s.*","",system(paste("wc -l", data_file), intern=TRUE)))
noquote(paste(data_file, "dimensions:", file_rows))
data <- read.table(data_file, sep="\t", header=TRUE, row.names=1, quote="\"*", comment.char = "")
# Check if data was loaded correctly
if (!  file_rows == nrow(data) + 1) {
  noquote(paste("Dataframe was not loaded correctly: number of rows in the dataframe (including header)", nrow(data)+1, "is not equal to number of rows in the file", file_rows))
  quit()
}

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
  
  
## Unify sample type for all non COVID-19
noquote("Unifying sample_type for all not COVID19_BALF, renaming to Another_BALF")
meta$sample_type[meta$sample_type != "COVID19_BALF"] <- "Another_BALF"

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

## COVID19_BALF vs Another_BALF
## Define contrasts, extract results table, and shrink the log2 fold changes
contrast_covid19 <- c("sample_type", "COVID19_BALF", "Another_BALF")
res_covid19_unshrunken <- results(dds, contrast=contrast_covid19, alpha = 0.05)
res_covid19 <- lfcShrink(dds, contrast=contrast_covid19, res=res_covid19_unshrunken, type="ashr")

## Save GO terms by padj values, add metadata for GO terms and write results
de_terms_covid19 <- res_covid19 %>% data.frame() %>% arrange(padj)
de_terms_covid19 <- merge(meta_go_term, de_terms_covid19, by="row.names") 
de_terms_covid19 <- arrange(de_terms_covid19, padj)
row.names(de_terms_covid19)<- de_terms_covid19$Row.names
de_terms_covid19 <- de_terms_covid19[-1]
go_term_meta_de_file = paste(output_dir, "de_go_terms_covid19_another.tsv", sep="")
noquote(paste("Saving differentially expressed GO terms in COVID19_BALF vs Another_BALF to", go_term_meta_de_file))
write.table(de_terms_covid19, file=go_term_meta_de_file, sep="\t", quote=F, col.names=NA)


## Select top 20 GO terms by padj values, add metadata for GO terms and write results
top20_terms_covid19 <- res_covid19 %>% data.frame() %>% arrange(padj) %>% head(n=20)
top20_terms_covid19 <- merge(meta_go_term, top20_terms_covid19, by="row.names") 
top20_terms_covid19 <- arrange(top20_terms_covid19, padj)
row.names(top20_terms_covid19)<- top20_terms_covid19$Row.names
top20_terms_covid19 <- top20_terms_covid19[-1]
go_term_meta_top20_file = paste(output_dir, "top20_go_terms_covid19_another.tsv", sep="")
noquote(paste("Saving top20 significant GO terms in COVID19_BALF vs Another_BALF to", go_term_meta_top20_file))
write.table(top20_terms_covid19, file=go_term_meta_top20_file, sep="\t", quote=F, col.names=NA)


# Plot expression for single term
noquote("Plotting top 20 significant GO terms")
for (term in rownames(top20_terms_covid19)) {
  plot_data <- plotCounts(dds, term, intgroup="sample_type", returnData = TRUE)
  plot_file = paste(output_dir, gsub(":", "-", term), ".png", sep="")
  png(file=plot_file, width=1200, height=900, res=120)
  plot_to_save <- ggplot(plot_data, aes(x=sample_type, y=count, color=sample_type)) + 
    ggtitle(paste(term, paste(meta_go_term[term,], collapse=", "), sep=", ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("normalized count") +
    geom_jitter(size=6) + 
    scale_y_log10()
  print(plot_to_save)
  dev.off()
}
noquote(paste("Plots for", paste(rownames(top20_terms_covid19), collapse=", "), "saved in", output_dir))

#Sys.sleep(120)

## Select top GO terms by depth, add metadata for GO terms and write results
top_depth_terms_covid19 <- res_covid19 %>% data.frame() %>% arrange(padj) %>% head(n=200)
top_depth_terms_covid19 <- merge(meta_go_term, top_depth_terms_covid19, by="row.names") 
top_depth_terms_covid19 <- top_depth_terms_covid19[top_depth_terms_covid19$depth >= 10, ]
row.names(top_depth_terms_covid19)<- top_depth_terms_covid19$Row.names
top_depth_terms_covid19 <- top_depth_terms_covid19[-1]
go_term_meta_top_depth_file = paste(output_dir, "top_depth_go_terms_covid19_another.tsv", sep="")
noquote(paste("Saving top significant GO terms with depth >=10 in COVID19_BALF vs Another_BALF to", go_term_meta_top_depth_file))
write.table(top_depth_terms_covid19, file=go_term_meta_top_depth_file, sep="\t", quote=F, col.names=NA)

# Plot terms with depth >=10 in top 200
noquote("Plotting significant GO terms with depth >=10 from top 200")
for (term in rownames(top_depth_terms_covid19)) {
  plot_data <- plotCounts(dds, term, intgroup="sample_type", returnData = TRUE)
  plot_file = paste(output_dir, gsub(":", "-", term), ".png", sep="")
  png(file=plot_file, width=1200, height=900, res=120)
  plot_to_save <- ggplot(plot_data, aes(x=sample_type, y=count, color=sample_type)) + 
    ggtitle(paste(term, paste(meta_go_term[term,], collapse="\n"), sep=", ")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("normalized count") +
    geom_jitter(size=6) + 
    scale_y_log10()
  print(plot_to_save)
  dev.off()
}
noquote(paste("Plots for", paste(rownames(top_depth_terms_covid19), collapse=", "), "saved in", output_dir))
