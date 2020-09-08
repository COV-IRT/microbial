library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)


# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
ps_balf_comp <- microbiome::transform(ps_balf, "compositional")

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(ps_balf_comp)
count <- as.matrix(t(dat))
fit <- lapply(1:10,dmn, count = count, verbose=TRUE)
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace
#plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

#pick the optimal model
best <- fit[[which.min(unlist(lplc))]]
#mixture parameters pi and theta
mixturewt(best)

#sample component assignments
ass <- apply(mixture(best), 1, which.max)
# Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}