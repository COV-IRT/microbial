library(microbiome)
library(DirichletMultinomial)
library(phyloseq)
#Set your
options(width=70, digits=2)
full <- FALSE
.qualitative <- DirichletMultinomial:::.qualitative
dev.off <- function(...) invisible(grDevices::dev.off(...))

#convert counts to a matrix
dat <- abundances(bac_pseq_prune)
count <- as.matrix(t(dat))



lvls <- c("Control_Healthy", "Control_Sick", "COVID19")
pheno<-factor(sample_data(bac_pseq_prune)$case, levels=lvls)
pheno

counts <- lapply(levels(pheno), csubset, count, pheno)
sapply(counts, dim)
#  [,1] [,2] [,3]
#  [1,]   32   61   33
#  [2,] 2077 2077 2077

keep <- c("Control_Healthy", "COVID19")
count <- count[pheno %in% keep,]
pheno <- factor(pheno[pheno %in% keep], levels=keep)

#The dmngroup function identifies the best (minimum Laplace score) 
#Dirichletmultinomial model for each group.
if (full) {
  bestgrp <- dmngroup(count, pheno, k=1:5, verbose=TRUE,
                        mc.preschedule=FALSE)
  save(bestgrp, file=file.path(tempdir(), "bestgrp.rda"))
  } else data(bestgrp)

'''The Control Health group is described by a model with one component, 
the Covid19 group by a model with three components. 
Three of the four Dirichlet components of the original single group 
(best) model are represented in the COVID19 group, the other in the Healthy group. 
The total Laplace score of the two group model is less than of the single-group model,
indicating information gain from considering groups separately.'''

bestgrp

lapply(bestgrp, mixturewt)



c(sapply(bestgrp, laplace),
  "Control_Healthy+COVID19"=sum(sapply(bestgrp, laplace)),
  Single=laplace(best))

"The predict function assigns samples to classes; the confusion matrix shows
that the classifier is moderately effective."

xtabs(~pheno + predict(bestgrp, count, assign=TRUE))

"The cvdmngroup function performs cross-validation. This is a computationally
expensive step"

if (full) 
{
     xval <- cvdmngroup(nrow(count), 
                        count, 
                        c(Control_Healthy=1, COVID19=3), 
                        pheno,
                        verbose=TRUE, 
                        mc.preschedule=FALSE)
     
     save(xval, file=file.path(tempdir(), "xval.rda"))
} 
else data(xval)

'''Figure 5 shows an ROC curve for the single and two-group classifier. 
The single group classifier is performing better than the two-group classifier.'''

bst <- roc(pheno[rownames(count)] == "COVID19",
              predict(bestgrp, count)[,"COVID19"])
bst$Label <- "Single"
two <- roc(pheno[rownames(xval)] == "COVID19",
              xval[,"COVID19"])
two$Label <- "Two group"
both <- rbind(bst, two)
pars <- list(superpose.line=list(col=.qualitative[1:2], lwd=2))
pdf("roc.pdf")
xyplot(TruePostive ~ FalsePositive, group=Label, both,
         + type="l", par.settings=pars,
         + auto.key=list(lines=TRUE, points=FALSE, x=.6, y=.1),
         + xlab="False Positive", ylab="True Positive")
dev.off()

#toLatex(sessionInfo())
"I. Holmes, K. Harris, and C. Quince. Dirichlet multinomial mixtures: Generative models for microbial metagenomics. PLoS ONE, 7(2):e30126, 02
2012."
