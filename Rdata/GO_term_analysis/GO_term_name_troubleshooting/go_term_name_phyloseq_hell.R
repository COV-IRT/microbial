"###################MIKE LEE PLEASE HELPPPP#############################
########WHY ARE THE GO_TERM NAMES BREAKING THE DATAFRAMES!?!?!#########"


#load libraries
library(phyloseq)

#make a little function that returns the top 1k 
#go term names from a phyloseq object

check_go_names<-function(x){
  test<-as.data.frame(tax_table(x))
  head(test$name, n=1000L)
}

#Ok Here is my original phyloseq object containing the bacterially derived
#raw biological processes counts and associated sample metadata
pseq<-readRDS("pseq.RDS")
pseq
check_go_names(pseq) # no problems here


# preprocessing
# this is fine
pseq_no_neg<-subset_samples(pseq, sample_type!="neg_control")
check_go_names(pseq_no_neg)# no problems here
pseq_no_neg<-subset_samples(pseq_no_neg, sample_type!="Unknown")
check_go_names(pseq_no_neg)# no problems here

#This seems to be breaking the go term names
pseq_prune <- prune_taxa(taxa_sums(pseq_no_neg) > 10, pseq)
#### WARNING, IF YOU DO THIS YOU WILL FREEZE THE CONSOLE AND MAKE IT LAGGY####
#check_go_names(pseq_prune)
pseq_prune <- prune_samples(sample_sums(pseq_prune) > 10, pseq_prune)
#check_go_names(pseq_prune)


#TOLD YA....