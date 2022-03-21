primerstats<-read.csv(file = "primer_stats.txt", header = T, sep = "\t", row.names = 1)
primerstats$SNV_count
lm(data = primerstats, formula = (SNV~SNP))



library(ggplot2);library(reshape2)
data<- melt(primerstats)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+xlab("Num. of SNPs")+ylab("Density")
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)+xlab("Num. of SNPs")+ylab("Num. of Primers")
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot(alpha=0.25)+xlab(NULL)+ylab("Num. of SNPs")


print(paste("min",min(primerstats$SNP_count), 
            "max",max(primerstats$SNP_count),
            "CI", mean_ci(primerstats$SNP_count, error.limit = "both"),sep = " "))
mean(primerstats$SNV_count)
which.max(primerstats$SNV_count)
primerstats[174,]
primerstats[222,]
boxplot(mean_ci(primerstats$SNV_count, error.limit = "both"))

library(tidygraph)

data<-
