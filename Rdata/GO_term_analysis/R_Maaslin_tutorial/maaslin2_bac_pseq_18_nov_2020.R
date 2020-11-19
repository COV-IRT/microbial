library(DESeq2)
sample_info_tab<-sample_data(bac_pseq_no_neg)
sample_info_tab_phy <- sample_data(sample_info_tab)
deseq_counts<-phyloseq_to_deseq2(physeq = bac_pseq_no_neg,design = ~ 1) 
deseq_counts_vst <- estimateSizeFactors(deseq_counts, type = "poscounts")
vst_trans_count_tab <- assay(deseq_counts_vst)

#YAAAAAAAAAAAAAAAASSSSSSSSSSSSS
vst_trans_count_tab2 <- limma::removeBatchEffect(vst_trans_count_tab, sample_info_tab$publication)
#IT FIXED THE BATCH EFFECT WORRRRRRRRRRRRRRRRRRRKED##############
as.tibble(vst_trans_count_tab2)
vst_count_phy <- otu_table(vst_trans_count_tab2, taxa_are_rows=T)
vst_tax_phy <- tax_table(bac_pseq_no_neg)
vst_physeq <- phyloseq(vst_count_phy, vst_tax_phy,sample_data(bac_pseq))

# 
# library(tidyverse)
# p<-psmelt(bac_pseq_prune)
# dmn_sum<-p%>%
#   select(Sample,OTU,dmn,case,Abundance)%>%
#   group_by(OTU,dmn,case)%>%
#   pivot_wider(id_cols = c(Sample,dmn,case), names_from = OTU,values_from = Abundance)
# dmn_sum
# write.table(dmn_sum, "dmn_sum.tsv",sep="\t")


#############MaAsLIN2#############
#BiocManager::install("Maaslin2")
#rm(list=ls())
#getwd()

#dir.create("R_Maaslin_tutorial") # Create a new directory
setwd("R_Maaslin_tutorial") # Change the current working directory 
getwd() #check if directory has been successfully changed

#Load MaAsLin2 package into the R environment
library(Maaslin2)



library(speedyseq)
df_input_data2<-data.frame(t(otu_table(vst_physeq)))
df_input_metadata2<-data.frame(sample_data(vst_physeq))

Maaslin2(
  input_data = df_input_data2,
  input_metadata = df_input_metadata2,
  output="covirt_bac_pseq_prune_maaslin2",
  min_abundance = 10,
  min_prevalence = 0.1,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  random_effects = c("publication"),
  fixed_effects = c("case"),
  correction = "BH",
  standardize = TRUE,
  cores = 8,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50
)

library(tidyverse)
met<-as_tibble(df_input_metadata2)
met
# input_data <- system.file(
#   'extdata','HMP2_taxonomy.tsv', package="Maaslin2")
# input_metadata <-system.file(
#   'extdata','HMP2_metadata.tsv', package="Maaslin2")
write.table(df_input_metadata2,"bac_pseq_prune_metadata.tsv",sep = "\t")
write.table(df_input_data2,"bac_pseq_prune_counts.tsv",sep = "\t")
fixed<-c('sample_type', 'sequenced.respiratory.pathogens	', 'case	', 'age	', 'sex	', 'temp_degC	', 'cough	', 'glucocorticoid_therapy	', 'abx_therapy	', 'antiviral_therapy	', 'Oxygen_therapy	', 'date_of_onset	', 'days_delayed_hospitalization	', 'SPO2.	', 'outcome	', 'days.after.onset	', 'origin	', 'hospitaladmit_date	', 'hospital.release.date	', 'viral_titre_perc	', 'CRP.mg.L..1.	', 'disease_severity	', 'BALF_SARS.CoV.2_Ct	', 'oral_swab_SARS.CoV.2_Ct	', 'anal_swab_SARS.CoV.2_Ct	', 'blood_ab_SARS.CoV.2	', 'blood_pcr_SARS.CoV.2	', 'cell_counts	', 'white.blood.count..WBC...109.L	', 'red.blood.count..RBC...1012.L	', 'Neutrophils...109.L.L	', 'lymphocyte...109.L	', 'monocyte...109.L	', 'haemoglobin.Hb..g.L	', 'Platelet...109.L	', 'Albumin..g.L	', 'aspartate.aminotransferase.AST..U.L	', 'alanine.aminotransferase.ALTU.L	', 'creatine.kinase.CK..U.L	', 'creatinine.kinase.MB.isoenzyme.CK.MB..U.L	', 'lactate.dehydrogenase.LDH..U.L	', 'UREA..mmol.L	', 'CREA...mol.L	', 'X.IL.6..pg.mL.	', 'IFN...pg.mL	', 'IL.10..pg.mL	', 'smoking_status	', 'smoking_pack_per_year	', 'confirmatory_test_done	', 'diagnosis_history	', 'collection_location	', 'reads')
random<-c( 'publication','bioproject	', 'sample_name	', 'isolate_name	', 'clinical.lab	', 'host	', 'sampling_site	', 'collection_date	', 'release_date	', 'collected_by	', 'submitted_by	', 'sequence_type	', 'accession	', 'dmn	')
colnames(df_input_metadata2))
list

#Interactions
df_input_metadata2$case<-as.character(df_input_metadata2$case)
df_input_metadata2$case_Control_Healthy<-(df_input_metadata2$case_modified == "Control_Healthy")*df_input_metadata2$case


df_input_metadata$UC_dysbiosis = (df_input_metadata$diagnosis_modified == "UC") *
  df_input_metadata$dysbiosis


fit_data <- Maaslin2(
  input_data = df_input_data2,
  input_metadata = df_input_metadata2,
  output = 'covirt_maaslin2', 
  transform = "LOG",
  # min_abundance = 10,
  # min_prevalence = 10,
  # normalization = 'TSS',
  correction = "BH",
  fixed_effects=c('sample_type','sequenced.respiratory.pathogens','case','age','sex','temp_degC','cough','glucocorticoid_therapy','abx_therapy','antiviral_therapy','Oxygen_therapy','date_of_onset','days_delayed_hospitalization','SPO2.','outcome','days.after.onset','origin','hospitaladmit_date','hospital.release.date','viral_titre_perc','CRP.mg.L..1.','disease_severity','BALF_SARS.CoV.2_Ct','oral_swab_SARS.CoV.2_Ct','anal_swab_SARS.CoV.2_Ct','blood_ab_SARS.CoV.2','blood_pcr_SARS.CoV.2','cell_counts','white.blood.count..WBC...109.L','red.blood.count..RBC...1012.L','Neutrophils...109.L.L','lymphocyte...109.L','monocyte...109.L','haemoglobin.Hb..g.L','Platelet...109.L','Albumin..g.L','aspartate.aminotransferase.AST..U.L','alanine.aminotransferase.ALTU.L','creatine.kinase.CK..U.L','creatinine.kinase.MB.isoenzyme.CK.MB..U.L','lactate.dehydrogenase.LDH..U.L','UREA..mmol.L','CREA...mol.L','X.IL.6..pg.mL.','IFN...pg.mL','IL.10..pg.mL','smoking_status','smoking_pack_per_year','confirmatory_test_done','diagnosis_history','collection_location','reads'),
  random_effects=c('publication','bioproject','sample_name','isolate_name','clinical.lab','host','sampling_site','collection_date','release_date','collected_by','submitted_by','sequence_type','accession','dmn'),
  standardize = T,
  cores = 8)


fit_data2 = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "demo_output2", 
  fixed_effects = c("diagnosis", "dysbiosis"),
  reference = c("diagnosis,nonIBD"))
