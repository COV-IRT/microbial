"###################MIKE LEE PLEASE HELPPPP#############################
########WHY ARE THE GO_TERM NAMES BREAKING THE DATAFRAMES!?!?!#########"

#load libraries
library(phyloseq)
library(tidyverse)

#OSF derived combined output from seqscreen
raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T)
dim(raw_df)
 # [1] 26966  2020

# this was supposed to have 47,234 lines (or 47,233 rows if read in with header):
  # wc -l Combined_BALF_GO_Terms.tsv
  #  47234 Combined_BALF_GO_Terms.tsv
# this is due to some having quotes in them that R is trying to interpret as field-quoting characters, e.g.:
  # problem will occur in the "name" column, where there can be quotes mixed in, e.g.:
     # GO:0006458	biological_process	3	'de novo' protein folding
       # or whatever the F this is
     # GO:1900549	biological_process	4	N',N'',N'''-triacetylfusarinine C metabolic process
       # even though both of those close properly, this one wouldn't
     # GO:0061146	biological_process	4	Peyer's patch morphogenesis
  # i found those peeking at the command line, e.g.: grep "'" Combined_BALF_GO_Terms.tsv | head
# that is likely causing one or more of these "name" entries to be thousands of rows long

# we can turn that off by setting the quoting characters to nothing, quote=""
raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T, quote="")
  # now we get an error that says
# Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
#   line 35939 did not have 2020 elements

   # telling us something is up with that row (or the one before or after it depending on how things are being counted)
  # peaking there at the command line:
  # sed -n '35938,35940p' Combined_BALF_GO_Terms.tsv | cut -f 1-4
    # GO:0000982	molecular_function	3	DNA-binding transcription factor activity, RNA polymerase II-specific
    # GO:0001133	molecular_function	3	DNA-binding transcription factor activity, RNA polymerase II-specific
    # GO:1904067	molecular_function	3	ascr#2 binding

  # ah, there is a '#' in that one, nice GO...
  # so it's reading that as a comment and then not able to finish that row, we can turn off commenting characters like so:
raw_df <- read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = "")
  # now the whole thing reads in with the expected amount of rows:
dim(raw_df)

  # it's super-valuable to check what we read in matches the expected number of rows and columns :)

# doing it same as you were from here on, just with those added arguments to the initial read.table call


raw<-as_tibble(read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T, quote = "", comment.char = ""))
raw
# A tibble: 47,233 x 2,020     # good so far now

#do a little regex and fix some stuff
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

#transform the raw table by type of count (euk, term, bac, arc)
df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","Bacteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value, values_fill=0)


####There are multiple processes and values for a single sample so you cant convert the sample to columns
#make individual tibbles for biological processes and molecular fxn
bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")

#make individual tibbles for each type (bac, euk, term, arc, vir, etc)
bio_bac<-bio%>%filter(type=="bac")%>%select(-type)
bio_term<-bio%>%filter(type=="term")%>%select(-type)
mol_bac<-mol%>%filter(type=="bac")%>%select(-type)
mol_term<-mol%>%filter(type=="term")%>%select(-type)

#subselect tibbles for only the counts and go terminology
bio_bac_counts<-bio_bac%>%select(-c(namespace,depth,name))
bio_bac_tax<-bio_bac%>%select(GO_term,namespace,depth,name)

#convert them to dataframes for downstream import to phylsoeq
bio_bac_counts<-data.frame(bio_bac_counts, row.names=1)
bio_bac_tax<-data.frame(bio_bac_tax, row.names=1)

bio_bac_counts_phy <- otu_table(bio_bac_counts, taxa_are_rows=TRUE)
bio_bac_tax_phy <- tax_table(as.matrix(bio_bac_tax), errorIfNULL=TRUE)

# checking nothing similar is happening with the sample info table (though probably not because it has quoting characters)
  # number rows (minus one if reading in as with a header)
# wc -l Combined_BALF_GO_Terms_metadata.txt
#      168 Combined_BALF_GO_Terms_metadata.txt
  # number of columns
# head -n 1 Combined_BALF_GO_Terms_metadata.txt | tr "\t" "\n" | wc -l
#        70

bio_bac_sam<-as.data.frame(read.table("Combined_BALF_GO_Terms_metadata.txt",header = T, sep = "\t",row.names = 1))
dim(bio_bac_sam)
  # [1] 167  70    # good

#a little regex to fix the stupid filename
rownames(bio_bac_sam)<-rownames(bio_bac_sam)%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")

# making physeq object

pseq <- phyloseq(bio_bac_counts_phy, bio_bac_tax_phy, sample_data(bio_bac_sam_phy))
##########



#make a little function that returns the top 1k
#go term names from a phyloseq object

# check_go_names<-function(x){
#   test<-as.data.frame(tax_table(x))
#   head(test$name, n=1000L)
# }

  ## the conversion to dataframe part of this was mkaing it take a long time:
start <- Sys.time() ; check_go_names(pseq) ; end <- Sys.time()
end - start # Time difference of 4.654899 mins

  # and i wasn't sure we needed it, so making one without it
check_go_names<-function(x){
    head(tax_table(x)[,3], n=1000L)
}

start <- Sys.time() ; check_go_names(pseq) ; end <- Sys.time()
end - start # Time difference of 0.03080297 secs

# preprocessing
# this is fine
pseq_no_neg<-subset_samples(pseq, sample_type!="neg_control")
check_go_names(pseq_no_neg)# no problems here
pseq_no_neg<-subset_samples(pseq_no_neg, sample_type!="Unknown")
check_go_names(pseq_no_neg)# no problems here

#This seems to be breaking the go term names

pseq_prune <- prune_taxa(taxa_sums(pseq_no_neg) > 10, pseq)
#### WARNING, IF YOU DO THIS YOU WILL FREEZE THE CONSOLE AND MAKE IT LAGGY####
check_go_names(pseq_prune)
pseq_prune <- prune_samples(sample_sums(pseq_prune) > 10, pseq_prune)
check_go_names(pseq_prune)

# doesn't seem to be a problem anymore :)

#TOLD YA.... <- haha, made me laugh

### additional note
  # i'm still not sure it's a good idea to have this table in here this way, i had thought at first you were doing our "domain" "phylum" "class" etc was like this:
    # e.g. for GO:0031647: https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0031647#term=history
    #                col1                   col2                        col3                                col4
    # GO:0019319     biological process     biological regulation       regulation of biological quality    regulation of protein stability

   ## which again, still wouldn't work perfectly, because not all of them are just linear like this one is, e.g. this one: https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0019319#term=history

  # but i don't think we want to use phyloseq tax functions on the current format of the tax table, with like the "depth" column in there,
  # i'm pretty sure that is telling phyloseq that everything under depth "1" for instance, is all under the same taxonomy, when depth "1" really is equivallent to the rank name (like "genus")

  # though i think we also might still want/need the count propagation up all parent terms, since you're able to spend time on figuring these things now,
  # i think maybe filtering based on the GO term counts (our OTU/ASV IDs here) rather than their associated taxonomy, either before making the phyloseq object,
  # or with filter_taxa() maybe? (you're sharper on phyloseq than me, despite being called filter_taxa(), i think this is operating on our units, not their taxonomy)
    # e.g

pseq
pseq_prune <- filter_taxa(pseq, function(x) sum(x) > 10, prune = TRUE)
  # that comes out slightly differently than what you were running above, i'm not sure why currently and can't really dig in so just sending you uncomplete input on this for you to deal with currently, ha

