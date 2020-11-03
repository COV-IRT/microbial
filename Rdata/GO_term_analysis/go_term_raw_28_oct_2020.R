library("tidyverse")
#OSF derived combined output from seqscreen
raw<-as_tibble(read.table("Combined_BALF_GO_Terms.tsv", sep = "\t", row.names = NULL, header = T))

raw
##OSF derived casey's datasheet
raw2<-as_tibble(read.table("combined_seqscreen_GO_summary_v2.tsv", sep = "\t", row.names = NULL, header = T))

# #filter out only the raw counts
# counts<-raw %>%
#   select(namespace,depth,name,ends_with("_counts"))%>%
#   gsub("_counts","", colnames(counts))
# counts
# term_counts<-raw %>%
#   select(namespace,depth,name,ends_with("_term_counts"))
# bac_counts<-raw %>%
#   select(namespace,depth,name,ends_with("_bac_counts"))
# 
# #Do a little regex in order to get the sample names to match the metadata
# colnames(counts)<-gsub("_counts","", colnames(counts))
# colnames(term_counts)<-gsub("_term_counts","",colnames(term_counts))
# colnames(bac_counts)<-gsub("_bac_counts","",colnames(bac_counts))
# 
# counts
# 
# df<-counts %>% 
#   unite("Total",ends_with("_term"), remove = FALSE)%>%
#   separate(Total, sep = "_",colnames(bac_counts), remove = F)%>%
#   unite("Archaea",ends_with("_arc"), remove = FALSE)%>%
#   separate(Archaea, sep = "_",colnames(bac_counts), remove = F)%>%
#   unite("Bacteria",ends_with("_bac"), remove = FALSE)%>%
#   separate(Bacteria, sep = "_",colnames(bac_counts), remove = F)%>%
#   unite("Eukarya",ends_with("_euk"), remove = FALSE)%>%
#   separate(Eukarya, sep = "_",colnames(bac_counts), remove = F)%>%
#   unite("Viridae",ends_with("_vir"), remove = FALSE)%>% 
#   separate(Viridae, sep = "_",colnames(bac_counts), remove = F)%>%
#   unite("Unclassified",ends_with("_NA_tax"), remove = FALSE)%>%
#   separate(Unclassified, sep = "_",colnames(bac_counts), remove = F)#%>%
# 
# counts
# 
# # See vignette("pivot") for examples and explanation
# 
# # Simplest case where column names are character data
# 
# 
# 

# # Slightly more complex case where columns have common prefix,
# # and missing missings are structural so should be dropped.
# counts<-raw %>%
#   select(namespace,depth,name,ends_with("_counts"))%>%
#   gsub("_counts","", colnames(counts))%>%
#   pivot_longer(
#     cols = ,ends_with("_term")
#     names_to = "Total",
#     names_prefix = "wk",
#     values_to = "rank",
#     values_drop_na = TRUE
#   )
# 
# 
colnames(raw)<-gsub("NA_tax","unclass", colnames(raw))%>%str_replace_all("NC1_SRR7796663", "NC1.SRR7796663")


df<-raw %>%
  select(GO_term,namespace,depth,name,ends_with("_counts"))%>%
  pivot_longer(cols = -c(GO_term,namespace,depth,name),
               names_to =  c("sample","type","abund"),#c("Total", "Archaea","Bacteria","Eukarya", "Viridae", "Unclassified"),
               names_pattern = "(.*)_(.*)_(.*)")%>%
  select(-abund)%>%
  filter(value>1)%>%
  pivot_wider(names_from = sample, values_from=value)
df
####There are multiple processes and values for a single sample so you cant convert the sample to columns

bio<-filter(df, namespace=="biological_process")
mol<-filter(df, namespace=="molecular_function")

bio_bac<-bio%>%filter(type=="bac")%>%select(-type)
bio_term<-bio%>%filter(type=="term")%>%select(-type)
mol_bac<-mol%>%filter(type=="bac")%>%select(-type)
mol_term<-mol%>%filter(type=="term")%>%select(-type)





