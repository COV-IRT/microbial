lvls <- c("Control_Healthy", "Control_Sick", "COVID19")
pheno<-factor(sample_data(bac_pseq_prune)$case, levels=lvls)
pheno


length(taxa_names(bac_pseq_prune))

#Fit the model

bestgrp <- dmngroup(count, pheno, k=1:5, verbose=TRUE,simplify = TRUE, mc.preschedule=FALSE)
bestgrp

c(sapply(bestgrp, laplace),
  "Control_Sick+COVID19"=sum(sapply(bestgrp, laplace)),
  Single=laplace(best))

xtabs(~pheno + predict(bestgrp, count, assign=TRUE))

xval<-cvdmngroup(nrow(count), 
                 count, 
                 c(Control_Healthy=1, Control_Sick=2, COVID19=3), 
                 pheno,
                 verbose=TRUE, 
                 mc.preschedule=FALSE)

save(xval, file=file.path(tempdir(), "xval.rda"))
xval
#fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)
fit<-dmngroup(dmn,count, k=1:5, simplify = TRUE,.lapply = parallel::mclapply)
#fit <- mclapply(1:8, dmn, count = count, verbose=TRUE)