

library(matrixTests)
case_sum<-psmelt(bac_pseq)%>%
  select(Sample,OTU,case,Abundance)%>%
  group_by(OTU,case)%>%
  pivot_wider(id_cols = c(Sample,case), names_from = OTU,values_from = Abundance)

dim(case_sum)
mat<-case_sum[,3:13848]
krus_case<-col_kruskalwallis(x =mat,g = case_sum$case)
welch_case<-col_oneway_welch(x=mat, g=case_sum$case)
krus_case_sig<-krus_case%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_case_sig<-welch_case%>%filter(pvalue<0.01)%>%arrange(pvalue)
dim(krus_case_sig)
dim(welch_case_sig)
write.table(krus_case,"krus_case.tsv",sep = "\t")
write.table(welch_case,"welch_case.tsv",sep = "\t")



#########################################
dmn_sum<-psmelt(bac_pseq)%>%
  select(Sample,OTU,dmn,Abundance)%>%
  group_by(OTU,dmn)%>%
  pivot_wider(id_cols = c(Sample,dmn), names_from = OTU,values_from = Abundance)


mat<-dmn_sum[,3:856]
krus_dmn<-col_kruskalwallis(x =mat,g = dmn_sum$dmn)
welc_dmnh<-col_oneway_welch(x=mat, g=dmn_sum$dmn)

krus_dmn_sig<-krus_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_dmn_sig<-welch_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_sig
dim(krus_dmn_sig)
dim(welch_dmn_sig)
krus_sig
write.table(krus,"krus_dmn.tsv",sep = "\t")

colnames(krus_sig)
colnames(krus_sig)
krus_sig_names<-intersect(rownames(krus_case_sig), rownames(krus_dmn_sig))
welch_sig_names<-intersect(rownames(welch_case_sig), rownames(welch_dmn_sig))
sig_names<-intersect(krus_sig_names,welch_sig_names)
#sig_names<-as_tibble(sig_names)
#colnames(sig_names)<-"names"




# library(matrixTests)
# 
# 
# case_sum<-psmelt(bac_pseq_prune_deeper2)%>%
#   select(Sample,OTU,case,Abundance)%>%
#   group_by(OTU,case)%>%
#   pivot_wider(id_cols = c(Sample,case), names_from = OTU,values_from = Abundance)
# 
# 
# mat<-case_sum[,3:856]
# krus_case<-col_kruskalwallis(x =mat,g = case_sum$case)
# welch_case<-col_oneway_welch(x=mat, g=case_sum$case)
# krus_case_sig<-krus_case%>%filter(pvalue<0.01)%>%arrange(pvalue)
# welch_case_sig<-welch_case%>%filter(pvalue<0.01)%>%arrange(pvalue)
# dim(krus_case_sig)
# dim(welch_case_sig)
# write.table(krus,"krus_case.tsv",sep = "\t")
# 
# 
# #########################################
# dmn_sum<-psmelt(bac_pseq_prune_deeper2)%>%
#   select(Sample,OTU,dmn,Abundance)%>%
#   group_by(OTU,dmn)%>%
#   pivot_wider(id_cols = c(Sample,dmn), names_from = OTU,values_from = Abundance)
# 
# 
# mat<-dmn_sum[,3:856]
# krus_dmn<-col_kruskalwallis(x =mat,g = dmn_sum$dmn)
# welc_dmnh<-col_oneway_welch(x=mat, g=dmn_sum$dmn)
# 
# krus_dmn_sig<-krus_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)
# welch_dmn_sig<-welch_dmn%>%filter(pvalue<0.01)%>%arrange(pvalue)
# welch_sig
# dim(krus_dmn_sig)
# dim(welch_dmn_sig)
# krus_sig
# write.table(krus,"krus_dmn.tsv",sep = "\t")
# 
# colnames(krus_sig)
# colnames(krus_sig)
# krus_sig_names<-intersect(rownames(krus_case_sig), rownames(krus_dmn_sig))
# welch_sig_names<-intersect(rownames(welch_case_sig), rownames(welch_dmn_sig))
# sig_names<-intersect(krus_sig_names,welch_sig_names)
# #sig_names<-as_tibble(sig_names)
# #colnames(sig_names)<-"names"
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# sig_bac_pseq_prune_deeper2<-prune_taxa(taxa = sig_names,x = bac_pseq_prune_deeper2)
# sig_bac_pseq_prune_deeper2
# count<-abundances(sig_bac_pseq_prune_deeper2)
# select <- order(rowMeans(count),decreasing=TRUE)
# select2<-sqrt((count)[select,])
# tmp<-rownames(select2)
# #select2$mol<-tmp
# select3<-as.table(select2)
# test<-t.test(select3)
# test
# library(matrixTests)
# library(genefilter)
# select3<-as_tibble(select2)%>%summarise(std=rowFtests(select2))%>% arrange(desc(std))
# select3
# select2
# select2$mol<-tmp
#
#
# select2
# select2<-as_tibble(select2)
# select2
# select2$mean<-rowMeans(select2)
# select3
# select2<-filter(.data = select2,rowMeans(select2)<10)
#
# sam<-data.frame(sample_data(sig_bac_pseq_prune_deeper2))
# df<-as.data.frame(sample_data(sig_bac_pseq_prune_deeper2))
# df<-as_tibble(df)
# df<-df%>%select(dmn, publication, sample_type,case)#,dmn,body_site)
# df
# df<-as.data.frame(df)
# df
#
# select2<-sqrt((count)[select,])
# dim(select2)
# colnames(select2) <- colnames(otu_table(bac_pseq_prune_deep))
# length(rownames(select2))
# length(row.names(df))
# row.names(df) <- colnames(select2)
#
#
# library(ggsci)
# mypal <- pal_aaas("default", alpha = 1)(10)
# mypal
# library("scales")
# library(RColorBrewer)
# library(viridis)
# library(pheatmap)
# df_row<-as.data.frame(fitted(best))
# df<-data.frame(df)
# df
# colnames(df)<-c("dmm_cluster", "Publication","Sample_Type","Case")
# # Specify colors
# ann_colors = list(
#   dmm_cluster=c("1"="#3B4992FF","2"="#EE0000FF","3"="#008B45FF","4"="#631879FF","5"="#008280FF"),
#   Publication=c("Xiong"="#008280FF", "Shen"="#631879FF","Michalovich"="#BB0021FF","Chen"="#3B4992FF","Wu"="#EE0000FF","Zhou"="orange","Ren"="#111111"),
#   Sample_Type =c("COVID_19"="firebrick","Healthy"="forestgreen","Community_acquired_pneumonia"="dodgerblue4",
#                  "Obese_Asthma"="goldenrod3",
#                  "Obese_Smoker"="goldenrod4",
#                  "Obese"="goldenrod1",
#                  "Asthma"= "darkorange2",
#                  "Asthma_Smoker"="darkorange4",
#                  "Asthma_Ex_smoker"="darkorange3",
#                  "Smoker"="gray27",
#                  "Obese_Asthma_Smoker"="black"),
#   Case=c("COVID19"="firebrick", "Control_Healthy"="forestgreen","Control_Sick"="dodgerblue4"))
#
#
# xx <- pheatmap(mat = select2,
#                color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(1000),
#                annotation_col=df,
#                annotation_colors = ann_colors,
#                clustering_distance_rows = "euclidean",
#                clustering_distance_cols = "euclidean",
#                annotation_row = df_row)
# dim(select2)


