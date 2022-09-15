library(base, quietly = TRUE)
library(methods, quietly = TRUE)
library(datasets, quietly = TRUE)
library(utils, quietly = TRUE)
library(grDevices, quietly = TRUE)
library(graphics, quietly = TRUE)
library(stats, quietly = TRUE)
library(stats4, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
library(S4Vectors, quietly = TRUE)
library(IRanges, quietly = TRUE)
library(DirichletMultinomial, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tibble, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(readr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(forcats, quietly = TRUE)
library(phyloseq, quietly = TRUE)
library(Maaslin2, quietly = TRUE)
library(microbiome, quietly = TRUE)
library(cartography, quietly = TRUE)
library(pheatmap, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(GenomeInfoDb, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(MatrixGenerics, quietly = TRUE)
library(Biobase, quietly = TRUE)
library(SummarizedExperiment, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(cluster, quietly = TRUE)
library(dendextend, quietly = TRUE)
library(rslurm, quietly = TRUE)


.rslurm_func <- readRDS('f.RDS')
.rslurm_params <- readRDS('params.RDS')
.rslurm_more_args <- readRDS('more_args.RDS')
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
.rslurm_istart <- .rslurm_id * 5 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 5, nrow(.rslurm_params))
.rslurm_result <- do.call(parallel::mcmapply, c(
    FUN = .rslurm_func,
    .rslurm_params[.rslurm_istart:.rslurm_iend, , drop = FALSE],
    MoreArgs = list(.rslurm_more_args),
    mc.cores = 2,
    mc.preschedule = TRUE,
    SIMPLIFY = FALSE))

saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
