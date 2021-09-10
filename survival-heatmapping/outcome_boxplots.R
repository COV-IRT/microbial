library(ggpubr)
counts<-as_tibble(samples_with_outcomes_vst_counts_df,rownames="term")
vec<-as_tibble(rownames(samples_with_outcomes_percentages_df))
samples_with_outcomes_percentages_df
vec
box_info<-as_tibble(GO_annots_heatmap_info,rownames="value")
box_info$value<-gsub("-containing"," containing",box_info$value)
vec<-full_join(vec,box_info)
vec$value2<-vec$value
#vec<-vec%>%separate(value,into = c("term","name"),sep = "-")

vec<-vec%>%separate(value2,into = c("term","name"),sep = "-")
vec
counts
samples_with_outcomes_percentages_df$value<-rownames(samples_with_outcomes_percentages_df)
samples_with_outcomes_percentages_df$value<-gsub("-containing"," containing",samples_with_outcomes_percentages_df$value)
samples_with_outcomes_percentages_df
df<-full_join(vec,counts)
vec
df<-full_join(vec,samples_with_outcomes_percentages_df)


dim(samples_with_outcomes_GO_for_heatmap)


a<-as_tibble(samples_with_outcomes_sample_annots_df,rownames = "ascession")
a


b<-df%>%pivot_longer(cols = -c(term,value,name,Depth,`Depth 1 Parent(s)`,NameSpace),names_to = "ascession",values_to = "vst_count")
b
b<-b[complete.cases(b),]
range(b$vst_count)# 1.438274 26.820147

d<-full_join(a,b)%>%filter(Outcome!="NA")%>%filter(!is.na(name))
d

unique(d$Outcome)
d$name<-as_factor(d$name)

my_comparisons<-list(c("Deceased","Survived"))
my_pal<-c("black","goldenrod1")
d$`Depth 1 Parent(s)`
ggviolin(data = d,
         x = "Outcome",
         y = "vst_count",
         add = "jitter",
         trim = T,
         color = "Outcome",
         orientation="vertical",
         palette = my_pal)+  
  stat_compare_means(size=3,comparisons = my_comparisons,method = "t.test")+
  facet_wrap(drop = T,scales = "free",facets= vars(`Depth 1 Parent(s)`,name),as.table = T,strip.position = "top") +
  theme(strip.background = element_blank(),text =  element_text(size=8),
        strip.placement = "outside")

t.test()

box_info<-as_tibble(GO_annots_heatmap_info)


################################################################################
cat<-full_join(a,b)%>%filter(Outcome!="NA")%>%filter(!is.na(name))
cat$Depth_Parent<-cat$`Depth 1 Parent(s)`
cat$Depth_Parent
cat<-cat%>%filter(Depth_Parent=="GO:0003824 | catalytic activity")%>%filter(!is.na(name))
cat<-cat%>%filter(name!="pyrophosphatase activity")
unique(d$Outcome)
cat$name<-as_factor(cat$name)
unique(cat$name)
my_comparisons<-list(c("Deceased","Survived"))
my_pal<-c("black","goldenrod1")

ggboxplot(data = cat,
         x = "Outcome",
         y = "vst_count",
         add = "jitter",
         trim = T,
         color = "Outcome",
         orientation="vertical",
         palette = my_pal)+  
  stat_compare_means(size=4,face="bold",comparisons = my_comparisons,method = "t.test")+
  facet_wrap(nrow = 1,drop = T,scales = "free",facets= ~name,as.table = T,strip.position = "top") +
  theme(strip.background = element_blank(),text =  element_text(size=12,face="bold"),
        strip.placement = "outside")

################################################################################

################################################################################
cat<-full_join(a,b)%>%filter(Outcome!="NA")%>%filter(!is.na(name))
cat$Depth_Parent<-cat$`Depth 1 Parent(s)`
cat<-cat%>%filter(!is.na(name))
unique(cat$name)
vec<-c("GO:0008152 | metabolic process","GO:0008152 metabolic process | GO:0009987  cellular process")
vec<-c("carbohydrate metabolic process" ,
       "organonitrogen compound catabolic process",
       "RNA metabolic process",
       "phosphorylation",
       "RNA phosphodiester bond hydrolysis",
       "nucleobase containing compound biosynthetic process")
cat<-cat%>%filter(name%in%vec)%>%filter(!is.na(name))
cat$Outcome
cat$name
unique(cat$name)

cat$name<-as_factor(cat$name)
my_comparisons<-list(c("Deceased","Survived"))
my_pal<-c("black","goldenrod1")

ggboxplot(data = cat,
          x = "Outcome",
          y = "vst_count",
          add = "jitter",
          trim = T,
          color = "Outcome",
          orientation="vertical",
          palette = my_pal)+  
  stat_compare_means(size=4,face="bold",comparisons = my_comparisons,method = "t.test")+
  facet_wrap(nrow = 1,drop = T,scales = "free",facets= ~name,as.table = T,strip.position = "top") +
  theme(strip.background = element_blank(),text =  element_text(size=11,face="bold"),
        strip.placement = "outside")

################################################################################
cat<-full_join(a,b)%>%filter(Outcome!="NA")%>%filter(!is.na(name))
cat$Depth_Parent<-cat$`Depth 1 Parent(s)`

cat<-cat%>%filter(!is.na(name))
unique(cat$name)
vec<-c("GO:0005488 | binding")
vec<-c("carbohydrate metabolic process" ,
       "organonitrogen compound catabolic process",
       "RNA metabolic process",
       "phosphorylation",
       "RNA phosphodiester bond hydrolysis",
       "nucleobase containing compound biosynthetic process")
cat<-cat%>%filter(Depth_Parent%in%vec)%>%filter(!is.na(name))
cat
#cat<-cat%>%filter(name%in%vec)%>%filter(!is.na(name))
cat$Outcome
cat$name
unique(cat$name)

cat$name<-as_factor(cat$name)
my_comparisons<-list(c("Deceased","Survived"))
my_pal<-c("black","goldenrod1")

ggboxplot(data = cat,
          x = "Outcome",
          y = "vst_count",
          add = "jitter",
          trim = T,
          color = "Outcome",
          orientation="vertical",
          palette = my_pal)+  
  stat_compare_means(size=4,face="bold",comparisons = my_comparisons,method = "t.test")+
  facet_wrap(nrow = 1,drop = T,scales = "free",facets= ~name,as.table = T,strip.position = "top") +
  theme(strip.background = element_blank(),text =  element_text(size=11,face="bold"),
        strip.placement = "outside")
