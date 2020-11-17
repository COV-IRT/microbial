library(matrixTests)


case_sum<-psmelt(mol_bac_pseq_prune_deep)%>%
  select(Sample,OTU,case,Abundance)%>%
  group_by(OTU,case)%>%
  pivot_wider(id_cols = c(Sample,case), names_from = OTU,values_from = Abundance)


case_sum

mat<-case_sum[,3:2387]
krus<-col_kruskalwallis(x =mat,g = case_sum$case)
welch<-col_oneway_welch(x=mat, g=case_sum$case)

welch
krus_sig<-krus%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_sig<-welch%>%filter(pvalue<0.01)%>%arrange(pvalue)
welch_sig
dim(krus_sig)
krus_sig
write.table(krus,"krus_case.tsv",sep = "\t")
krus_sig

