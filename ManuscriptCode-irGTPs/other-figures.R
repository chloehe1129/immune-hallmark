############### Study figure 4 cluster 2 ###############
library(dplyr)
library(tidyr)
library(tibble)
library(networkD3)
#### River Plot ###
meta$Immune.Subtype_archetype = as.factor(kmean_cluster_archetype$cluster)
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
data_long = as.data.frame(table(meta$Immune.Subtype_new,meta$Immune.Subtype_archetype))
colnames(data_long) <- c("source", "target", "value")
# data_long$source = ifelse(data_long$source=="1","Immune1",
#                           ifelse(data_long$source=="2","Immune2",
#                                  ifelse(data_long$source=="3","Immune3","Immune4")))
data_long$source = ifelse(data_long$source=="1","Immune1",
                          ifelse(data_long$source=="2","Immune2",
                                 ifelse(data_long$source=="3","Immune3",
                                        ifelse(data_long$source=="4","Immune4",
                                               ifelse(data_long$source=="5","Immune5","Immune6")))))
data_long$target = ifelse(data_long$target=="1","Archetype1",
                          ifelse(data_long$target=="2","Archetype2",
                                 ifelse(data_long$target=="3","Archetype3",
                                        ifelse(data_long$target=="4","Archetype4",
                                               ifelse(data_long$target=="5","Archetype5","Archetype6")))))
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1
sankeyNetwork(Links = data_long, Nodes = nodes,Source = "IDsource",
              Target = "IDtarget",Value = "value", NodeID = "name", height=300,width=350,
              sinksRight=T,nodeWidth=20,fontSize=15,nodePadding=5,LinkGroup = 'source')


#### River Plot ###
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
meta_twoless=meta[-which(meta$Immune.Subtype%in%c("C4","C6")),]
data_long = as.data.frame(table(meta_twoless$Immune.Subtype,meta_twoless$Immune.Subtype_new))
colnames(data_long) <- c("source", "target", "value")
data_long$source = ifelse(data_long$source=="C1","C1:Wound.Healing",
                          ifelse(data_long$source=="C2","C2:IFN-gamma",
                                 ifelse(data_long$source=="C3","C3:Inflammatory","C5:Immunologically.Quiet")))
# data_long$source = ifelse(data_long$source=="C1","C1:Wound.Healing",
#                           ifelse(data_long$source=="C2","C2:IFN-gamma",
#                                  ifelse(data_long$source=="C3","C3:Inflammatory",
#                                         ifelse(data_long$source=="C4","C4:Lymphocyte.Depleted",
#                                                ifelse(data_long$source=="C5","C5:Immunologically.Quiet","TGF.beta")))))
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1
sankeyNetwork(Links = data_long, Nodes = nodes,Source = "IDsource",
              Target = "IDtarget",Value = "value", NodeID = "name", height=300,width=350,
              sinksRight=T,nodeWidth=20,fontSize=15,nodePadding=5,LinkGroup = 'source')

gene_of_interest = c("IFNG","JAK1","STAT1","S100A1","E2F1","APC")
gene_df = as.data.frame(t(count[which(rownames(count)%in%gene_of_interest),]))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
gene_df = as.data.frame(apply(gene_df,2,range01))
gene_df = gene_df[rownames(meta),]
gene_df$Immune.Subtype_new = meta$Immune.Subtype_new
gene_df_long = melt(gene_df,id.var = "Immune.Subtype_new")
gene_df_long = gene_df_long[which(gene_df_long$Immune.Subtype_new%in%c("1","3","6")),]
gene_df_long$variable = factor(gene_df_long$variable,levels=c("IFNG","JAK1","STAT1","S100A1","E2F1","APC"))
gene_df_long$value = scale(gene_df_long$value,center = F,scale=T)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_cluster3_transition.pdf",
    width=5,height=4)
ggplot(gene_df_long,aes(x = Immune.Subtype_new,y=value,fill=Immune.Subtype_new))+
  geom_violin()+facet_wrap(~variable,scale="free",ncol=3)+
  theme(text = element_text(size = 30))+ylab("Expression")+theme_classic()+
  theme(legend.position = "none")+ stat_summary(fun.data="mean_sdl",geom="crossbar", width=0.2)+
  stat_compare_means(comparisons = list(c("1","3"),c("1","6"),c("3","6")),label = "p.signif")+xlab("Kmean Clusters")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
dev.off()

gene_of_interest = c("CD274","CTLA4","TIGIT")
gene_df = as.data.frame(t(count[which(rownames(count)%in%gene_of_interest),]))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
gene_df = as.data.frame(apply(gene_df,2,range01))
gene_df = gene_df[rownames(meta),]
gene_df$Immune.Subtype_new = meta$Immune.Subtype_new
gene_df_long = melt(gene_df,id.var = "Immune.Subtype_new")
gene_df_long = gene_df_long[which(gene_df_long$Immune.Subtype_new%in%c("1","3","6")),]
gene_df_long$value = scale(gene_df_long$value,center = T,scale=T)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_cluster136_checkpoint.pdf",
    width=5,height=2)
ggplot(gene_df_long,aes(x = Immune.Subtype_new,y=value,fill=Immune.Subtype_new))+
  geom_violin()+facet_wrap(~variable,scale="free",ncol=6)+
  theme(text = element_text(size = 30))+ylab("Expression")+theme_classic()+
  theme(legend.position = "none")+ stat_summary(fun.data="mean_sdl",size = 0.05)+
  stat_compare_means(comparisons = list(c("1","3"),c("1","6"),c("3","6")),label = "p.signif")+
  xlab("Kmean Clusters")+ylim(c(-1,3))+ annotate(geom="text", x=1.1, y=2, label="***",
                                                 color="black",hjust=0)
dev.off()

meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
df = as.data.frame(t(table(meta$Immune.Subtype,meta$Immune.Subtype_new)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("TCGA.Subtype","Immune.Subtype.new","Proportion")
n <- 4
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_cluster2_cancertype.pdf",width=9,height=6)
ggplot(prop, aes(fill=TCGA.Subtype, y=Proportion, x=Immune.Subtype.new)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("TCGA subtypes for each kmeans cluster")+
  scale_fill_manual(values=sample(col_vector, 6))+theme_classic()+xlab("Kmean Clusters")
dev.off()

meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
df = as.data.frame(t(table(meta$TCGA.Study,meta$Immune.Subtype_new)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("Cancer.Type","Immune.Subtype.new","Proportion")
n <- 33
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_cluster2_cancertype.pdf",width=9,height=6)
ggplot(prop, aes(fill=Cancer.Type, y=Proportion, x=Immune.Subtype.new)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Cancer composition for each kmeans cluster")+
  scale_fill_manual(values=sample(col_vector, 33))+theme_classic()+xlab("Kmean Clusters")
dev.off()

gene_of_interest = c("GZMK","CD8A","PRF1","CTLA4","PDCD1","TIGIT")
gene_df = as.data.frame(t(count[which(rownames(count)%in%gene_of_interest),]))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
gene_df = as.data.frame(apply(gene_df,2,range01))
gene_df = gene_df[rownames(meta),]
gene_df$Immune.Subtype_new = meta$Immune.Subtype_new
gene_df_long = melt(gene_df,id.var = "Immune.Subtype_new")
gene_df_long$variable = factor(gene_df_long$variable,levels=c("GZMK","CD8A","PRF1","CTLA4","PDCD1","TIGIT"))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_cluster2_exh_cyto.pdf",
    width=6,height=3)
ggplot(gene_df_long,aes(x = Immune.Subtype_new,y=value,fill=Immune.Subtype_new))+
  geom_violin(trim=T)+facet_wrap(~variable,scale="free",ncol=3)+xlab("Kmean Clusters")+
  coord_cartesian(ylim=c(0,0.05))+stat_compare_means(label='p.signif')+
  theme(text = element_text(size = 15))+ylab("Expression")+theme_classic()+
  theme(legend.position = "none")+ 
  annotate(geom="text", x=1.35, y=0.04, label="p-value<0.05",color="black",hjust=0)
dev.off()

ct_fraction = read.delim("/Volumes/she4/hallmark/data/validation-data/TCGA/TCGA.Kallisto.fullIDs.cibersort.relative.tsv")
ct_fraction = ct_fraction[!duplicated(ct_fraction$SampleID),]
ct_fraction$tcga_sample = substr(ct_fraction$SampleID,1,12)
ct_fraction = ct_fraction[!duplicated(ct_fraction$tcga_sample),]
rownames(ct_fraction) = gsub("\\.","-",ct_fraction$tcga_sample)
ct_fraction = ct_fraction[,c(3:24)]
# normalize <- function(x, na.rm = TRUE) x*100/sum(x)
# ct_fraction = as.data.frame(t(apply(ct_fraction,1,normalize)))
ct_fraction$sampleID = rownames(ct_fraction)
meta_fraction = merge(ct_fraction,meta[,c("Immune.Subtype_new","sampleID")],by="sampleID")
meta_fraction = meta_fraction[,-1]
tmp = aggregate(meta_fraction,by=list(meta_fraction$Immune.Subtype_new),FUN=mean)
tmp = tmp[,-c(1,ncol(tmp))]
df = data.frame(group = rep(colnames(tmp),6),
                cluster = rep(c(1:6),each=22),
                value = c(as.vector(unlist(tmp[1,])),
                          as.vector(unlist(tmp[2,])),
                          as.vector(unlist(tmp[3,])),
                          as.vector(unlist(tmp[4,])),
                          as.vector(unlist(tmp[5,])),
                          as.vector(unlist(tmp[6,]))))
ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+facet_wrap(~cluster)

ssGSEA_MP = read.csv(paste0(path_result,"TCGA_ssGSEA.csv"),row.names=1)
ssGSEA_MP = ssGSEA_MP[rownames(meta),]
meta = as.data.frame(cbind(ssGSEA_MP,meta))
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)

## ANOVA - Hallmark
# meta_hallmark = as.data.frame(cbind(meta[,58:107],meta$Immune.Subtype_new))
# p = NULL
# for(h in 1:50){
#   p = c(p,summary(aov(meta_hallmark[,h]~`meta$Immune.Subtype_new`,data = meta_hallmark))[[1]][1,5])
# }
Hallmark = read.csv("/Volumes/she4/hallmark/data/hallmark.csv")
temp = aggregate(meta,by=list(meta$Immune.Subtype_new),FUN = mean,na.rm=T)
temp = temp[,c(59:108)]
temp = temp[,Hallmark[,1]]
rownames(temp) = temp$Group.1
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4","5","6")
ann_row <- as.data.frame(read.csv("/Volumes/she4/hallmark/data/hallmark.csv")[,2])
names(ann_row) = "Hallmark.Group"
ann_row$Hallmark.Group = factor(ann_row$Hallmark.Group)
colnames(ann_row) <- c("Hallmark.Group")
set.seed(133)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
row = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))),8)
names(row) = unique(ann_row$Hallmark.Group)
row = list(row)
names(row) =c("Hallmark.Group")
set.seed(133)
rowAnn <- rowAnnotation(#df = ann_row,
  #annotation_width = unit(c(1, 4), 'cm'),
  foo = anno_block(gp = gpar(fill = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))),8))),
  gap = unit(1, 'mm'),show_annotation_name=F,labels = ann_row$Hallmark.Group)
h=Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = F,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Hallmark expression across immune subtypes",
        row_names_max_width = unit(12, "cm"),
        left_annotation = rowAnn,
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        })
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_subtypeFIG4_Hallmark.pdf",
    width=10,height=13)
draw(h, heatmap_legend_side="right")
dev.off()

meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
temp = aggregate(meta,by=list(meta$Immune.Subtype_new),FUN = mean,na.rm=T)
temp = temp[,c(31:58)]
rownames(temp) = temp$Group.1
temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4","5","6")
pdf("~/desktop/TCGA_novel_subtypeFIG4_KEGG.pdf",width=10,height=10)
Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "KEGG expression across different immune subtypes",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

## ANOVA - CellType-Signature
# meta_sig = as.data.frame(cbind(meta[,169:202],meta$Immune.Subtype_new))
# p = NULL
# for(h in 1:34){
#   p = c(p,summary(aov(meta_sig[,h]~`meta$Immune.Subtype_new`,data = meta_sig))[[1]][1,5])
# }
temp = aggregate(meta,by=list(meta$Immune.Subtype_new),FUN = mean,na.rm=T)
temp = temp[,c(170:205)]
temp = temp[,c(1:3,23:29,8:9,22,10:11,12,16:21,14:15)]
rownames(temp) = temp$Group.1
#temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4","5","6")
h=Heatmap(as.matrix(t(temp)), #[c(10,34,31,19),]
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Celltype across immune subtypes",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        },heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(6, "cm")))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_novel_subtypeFIG4_immunesig1_ordered.pdf",
    width=5.2,height=6)
draw(h, heatmap_legend_side="bottom")
dev.off()

temp = aggregate(meta,by=list(meta$Immune.Subtype_new),FUN = mean,na.rm=T)
temp = temp[,c(145:168)]
rownames(temp) = temp$Group.1
temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4","5","6")
pdf("~/desktop/TCGA_novel_subtypeFIG4_immunesig1.pdf",width=10,height=10)
Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "TCGA Immune expression across different immune subtypes",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

############### TCR-anchoring MP ###############

library(ggplot2)
enrich_L = readRDS(paste0("/Volumes/she4/hallmark/results/","c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_enrich_result.rds"))
tmp = enrich_L[[2]]
tmp = tmp[order(-tmp$Count,tmp$p.adjust),]
tmp = tmp[1:8,]
tmp$Description = factor(tmp$Description,levels=rev(tmp$Description))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/tcr_anchor.pdf",
    width=13,height=7)
ggplot(tmp,aes(x = Description,y=Count)) + 
  geom_bar(aes(fill = p.adjust),stat = "identity",position = "dodge") + 
  scale_y_log10()+coord_flip()+scale_fill_gradient(low="rosybrown2",high="paleturquoise3")+
  theme(text = element_text(size = 20))  + xlab("Pathways") + ylab("Core Enrichment")+
  theme(text = element_text(face = "bold")) + ggtitle("L_MP2: TCR Anchoring ORA")
dev.off()


library(ggpubr)
tmp = as.data.frame(cbind(meta[,"Lymphocyte.Infiltration.Signature.Score"],meta$Immune.Subtype))
tmp = tmp[which(tmp$V2%in%c("C1","C2","C3","C5")),]
tmp$V1 = as.numeric(tmp$V1)
tmp$V2 = as.factor(tmp$V2)
names(tmp)[2] = "Subtype"
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_C3-5_violin_all.pdf",
    width=4,height=2,onefile = F)
ggplot(tmp,aes(x=Subtype,y=V1,fill=Subtype))+geom_violin()+
  stat_compare_means(method="t.test",comparisons = list(c("C3","C5")),label="p.signif")+
  ylab("Lymphocyte Infiltration")+xlab("TCGA Subtypes")+
  scale_fill_manual(values=c("brown2","purple","cornflowerblue","#E69F00"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
dev.off()

library(ggpubr)
tmp = as.data.frame(cbind(meta[,151],kmean_cluster$cluster))
tmp = tmp[which(tmp$V2%in%c("2","4","5")),]
tmp$V1 = as.numeric(tmp$V1)
tmp$V2 = as.factor(tmp$V2)
names(tmp)[2] = "Subtype"
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_C3-5_violin.pdf",width=4,height=3,onefile = F)
ggplot(tmp,aes(x=Subtype,y=V1,fill=Subtype))+geom_violin()+
  stat_compare_means(method="wilcox.test",comparisons = list(c("2","4"),c("2","5"),c("4","5")))+
  ylab("Lymphocyte Infiltration")+xlab("TCGA Subtypes")+
  scale_fill_manual(values=c("cornflowerblue","#E69F00","green"))
dev.off()

############### Study figure 6 KM Curves ###############

library(survminer)
set.seed(1129)
kmean_cluster = kmeans(ssGSEA_MP, 6)
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
meta1 = meta[-which(is.na(meta$OS.time)),]
# meta1 = meta1[-which(meta1$OS.time==0),]
meta1$Immune_Subtypes = ifelse(meta1$Immune.Subtype=="C1","C1:Wound Healing",
                             ifelse(meta1$Immune.Subtype=="C2","C2:IFN-gamma",
                                    ifelse(meta1$Immune.Subtype=="C3","C3:Inflammatory",
                                           ifelse(meta1$Immune.Subtype=="C4","C4:Lymphocyte Depleted",
                                                  ifelse(meta1$Immune.Subtype=="C5","C5:Immunologically Quiet","C6:TGF-beta")))))
meta1 = meta1[which(meta1$Immune_Subtypes%in%c("C1:Wound Healing","C2:IFN-gamma","C3:Inflammatory","C5:Immunologically Quiet")),]
fit<- survfit(Surv(meta1$OS.time, meta1$OS.x) ~ Immune.Subtype, data = meta1)
a = ggsurvplot(fit, data = meta1,palette=c("brown2","purple","cornflowerblue","#E69F00"))+ggtitle("Kaplan Meier for TCGA Pan-cancer")+
  guides(colour = guide_legend(nrow = 1))
meta2 = meta1[which(meta1$Immune.Subtype %in% c("C3","C5")),]
meta3 = meta1[which(meta1$Immune.Subtype %in% c("C1","C2")),]
meta23 = meta1[which(meta1$Immune.Subtype %in% c("C1","C2","C3","C5")),]
res=pairwise_survdiff(Surv(meta2$OS.time, meta2$OS.x) ~ Immune.Subtype, data = meta2)
res2=pairwise_survdiff(Surv(meta3$OS.time, meta3$OS.x) ~ Immune.Subtype, data = meta3)
meta23$group = ifelse(meta23$Immune.Subtype%in%c("C1","C2"),"group1","group2")
res23=pairwise_survdiff(Surv(meta23$OS.time, meta23$OS.x) ~ group, data = meta23)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_KM_plot_subtypes.pdf",
    width=7,height=5,onefile = F)
a$plot + annotate("text",x = 3500, y = 0.8, size = 5,label = paste0("C3 vs C5: p value=",round(res$p.value,3)))
dev.off()
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_KM_plot_subtypes_pvalue.pdf",
    width=8,height=5,onefile = F)
a$plot + annotate("text",x = 1600, y = 0.2, size = 5,label = paste0("C3 vs C5: p value=",round(res$p.value,3)))+
  annotate("text",x = 1600, y = 0.1, size = 5,label = paste0("C1 vs C2: p value=",round(res2$p.value,3)))+
  annotate("text",x = 1800, y = 0, size = 5,label = paste0("C1+C2 vs C3+C5: p value=",round(res23$p.value,3)))
dev.off()

library(survminer)
set.seed(1129)
kmean_cluster = kmeans(ssGSEA_MP, 6)
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
meta1 = meta[-which(is.na(meta$OS.time)),]
fit<- survfit(Surv(meta1$OS.time, meta1$OS.x) ~ Immune.Subtype_new, data = meta1)
a = ggsurvplot(fit, data = meta1,conf.int=T,risk.table = T,
               risk.table.y.text.col = T, risk.table.y.text = FALSE)+
  ggtitle("Kaplan Meier across kmean clusters")+
  guides(colour = guide_legend(nrow = 2))
meta0 = meta1[which(meta1$Immune.Subtype_new %in% c("4","5")),]
meta2 = meta1[which(meta1$Immune.Subtype_new %in% c("2","5")),]
meta3 = meta1[which(meta1$Immune.Subtype_new %in% c("2","4")),]
meta4 = meta1[which(meta1$Immune.Subtype_new %in% c("1","3")),]
meta5 = meta1[which(meta1$Immune.Subtype_new %in% c("1","6")),]
meta6 = meta1[which(meta1$Immune.Subtype_new %in% c("3","6")),]
meta_groups = meta1
meta_groups$groups = ifelse(meta1$Immune.Subtype_new%in%c("1","3","6"),"Group1","Group2")
res0=pairwise_survdiff(Surv(meta0$OS.time, meta0$OS.x) ~ Immune.Subtype_new, data = meta0)
res=pairwise_survdiff(Surv(meta2$OS.time, meta2$OS.x) ~ Immune.Subtype_new, data = meta2)
res2=pairwise_survdiff(Surv(meta3$OS.time, meta3$OS.x) ~ Immune.Subtype_new, data = meta3)
res3=pairwise_survdiff(Surv(meta4$OS.time, meta4$OS.x) ~ Immune.Subtype_new, data = meta4)
res4=pairwise_survdiff(Surv(meta5$OS.time, meta5$OS.x) ~ Immune.Subtype_new, data = meta5)
res5=pairwise_survdiff(Surv(meta6$OS.time, meta6$OS.x) ~ Immune.Subtype_new, data = meta6)
res_groups=pairwise_survdiff(Surv(meta_groups$OS.time, meta_groups$OS.x) ~ groups, data = meta_groups)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_KM_plot_novelsubtypes.pdf",
    width=7,height=5,onefile = F)
# pdf("~/desktop/TCGA_KM_plot_novelsubtypes.pdf",
#     width=7,height=5,onefile = F)
a$plot + annotate("text",x = 8000, y = 1.1, size = 4,label = paste0("*2 vs 5: p value = ",formatC(res$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 1.0, size = 4,label = paste0("*2 vs 4: p value = ",formatC(res2$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.9, size = 4,label = paste0("*4 vs 5: p value = ",formatC(res0$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.8, size = 4,label = paste0(" 1 vs 3: p value = ",formatC(res3$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.7, size = 4,label = paste0(" 1 vs 6: p value = ",formatC(res4$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.6, size = 4,label = paste0("*3 vs 6: p value = ",formatC(res5$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 2500, y = 0.1, size = 4,label = paste0("*245 vs 136: p value = ",formatC(res_groups$p.value, format = "e", digits = 2)  ))
dev.off()

meta1$Immune.Subtype_new = factor(meta1$Immune.Subtype_new,levels = c("2","5","4","6","1","3"))
fit<- survfit(Surv(meta1$OS.time, meta1$OS.x) ~ Immune.Subtype_new, data = meta1)
a = ggsurvplot(fit, data = meta1,conf.int=T,surv.median.line="v",risk.table = T,
               risk.table.y.text.col = T, risk.table.y.text = FALSE,
               palette = c("goldenrod","dodgerblue","darkturquoise","magenta","salmon","green3"))+
  ggtitle("Kaplan Meier across kmean clusters")+
  guides(colour = guide_legend(nrow = 2))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_KM_plot_novelsubtypes_table.pdf",
    width=8,height=2,onefile = F)
a$table+ggtitle("k-mean clusters ordered from best to worst survival")
dev.off()

library(survminer)
set.seed(1129)
kmean_cluster = kmeans(ssGSEA_MP, 6)
meta$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
meta1 = meta[-which(is.na(meta$OS.time)),]
meta1_136 = meta1[which(meta1$Immune.Subtype_new%in%c("1","3","6")),]
fit<- survfit(Surv(meta1_136$OS.time, meta1_136$OS.x) ~ Immune.Subtype_new, data = meta1_136)
a = ggsurvplot(fit, data = meta1_136,conf.int=T,risk.table = T,
               risk.table.y.text.col = T, risk.table.y.text = FALSE)+
  ggtitle("Kaplan Meier across kmean clusters")+
  guides(colour = guide_legend(nrow = 1))
meta4_136 = meta1_136[which(meta1_136$Immune.Subtype_new %in% c("1","3")),]
meta5_136 = meta1_136[which(meta1_136$Immune.Subtype_new %in% c("1","6")),]
meta6_136 = meta1_136[which(meta1_136$Immune.Subtype_new %in% c("3","6")),]
res3_136=pairwise_survdiff(Surv(meta4_136$OS.time, meta4_136$OS.x) ~ Immune.Subtype_new, data = meta4_136)
res4_136=pairwise_survdiff(Surv(meta5_136$OS.time, meta5_136$OS.x) ~ Immune.Subtype_new, data = meta5_136)
res5_136=pairwise_survdiff(Surv(meta6_136$OS.time, meta6_136$OS.x) ~ Immune.Subtype_new, data = meta6_136)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_KM_plot_novelsubtypes136.pdf",
    width=9,height=4,onefile = F)
# pdf("~/desktop/TCGA_KM_plot_novelsubtypes.pdf",
#     width=7,height=5,onefile = F)
a$plot + annotate("text",x = 8000, y = 0.8, size = 4,label = paste0(" 1 vs 3: p value = ",formatC(res3_136$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.6, size = 4,label = paste0(" 1 vs 6: p value = ",formatC(res4_136$p.value, format = "e", digits = 2)  ))+
  annotate("text",x = 8000, y = 0.7, size = 4,label = paste0("*3 vs 6: p value = ",formatC(res5_136$p.value, format = "e", digits = 2)  ))
dev.off()

dist = as.data.frame(t(table(meta1$Immune_Subtypes,meta1$Immune.Subtype_new)))
dist = reshape2::acast(dist, Var1~Var2, value.var="Freq")
prop = apply(dist,2,prop.table)
prop = melt(prop)
names(prop) = c("Kmean_Cluster","TCGA_subtypes","Proportion")
prop$Kmean_Cluster = as.factor(prop$Kmean_Cluster)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_subtypes_in_KmeanClusters.pdf",width=6,height=4)
ggplot(prop, aes(fill=TCGA_subtypes, y=Proportion, x=Kmean_Cluster)) + xlab("Kmean Clusters") +
  geom_bar(position="fill", stat="identity")+ggtitle("TCGA subtype abundance for each kmean cluster")+
  theme(axis.text.x = element_text(face="bold", color=c("orchid3", "gold3", "chartreuse4","lightseagreen","skyblue3","salmon")))+
  scale_fill_manual(values=c("brown2","purple","cornflowerblue","#E69F00"))+theme_classic()
dev.off()

dat = as.data.frame(cbind(tcga.umap$layout,meta$type))
dat$V1 = as.numeric(dat$V1)
dat$V2 = as.numeric(dat$V2)
dat$V3 = as.factor(dat$V3)
dat = dat[-which(dat$V3=="THYM"),]
ggplot(dat,aes(x=V1,y=V2,colour=V3))+geom_point()+xlab("UMAP 1")+ylab("UMAP 2")

dist = as.data.frame(t(table(meta$Immune.Subtype,meta$Immune.Subtype_new)))
dist = reshape2::acast(dist, Var1~Var2, value.var="Freq")
prop = apply(dist,2,prop.table)
prop = melt(prop)
names(prop) = c("Cluster","TCGA_subtypes","Proportion")
prop$Cluster = as.factor(prop$Cluster)
pdf("~/desktop/stacked_celltype.pdf",width=8,height=6)
ggplot(prop, aes(fill=TCGA_subtypes, y=Proportion, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")+ggtitle("TCGA subtypes Abundance for each cluster")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

############### Study figure 6E - PI clusters ###############
test_meta_subset = test_meta[-which(test_meta$clinical_stage%in%c("[Not Available]",
                                                                  "[Discrepancy]",
                                                                  "[Not Applicable]","I","IIa","IIb","III","IVa","IVb")),]
test_meta_subset$Clinical.Stage = ifelse(test_meta_subset$clinical_stage%in%c("Stage I","Stage IA","Stage IA1",
                                                                              "Stage IB","Stage IB1","Stage IB2",
                                                                              "Stage IC"),"Stage I",
                                         ifelse(test_meta_subset$clinical_stage%in%c("Stage II","Stage IIA",
                                                                                     "Stage IIA1","Stage IIA2",
                                                                                     "Stage IIB","Stage IIC"),"Stage II",
                                                ifelse(test_meta_subset$clinical_stage%in%c("Stage III","Stage IIIA",
                                                                                            "Stage IIIB","Stage IIIC",
                                                                                            "Stage IIIC1","Stage IIIC2",
                                                                                            "Stage IS"),"Stage III","Stage IV")))
df = as.data.frame(t(table(test_meta_subset$Clinical.Stage,test_meta_subset$risk_group)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("Clinical.Stage","Risk.Group","Proportion")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_clinicalstage.pdf",
    width=4,height=3)
ggplot(prop, aes(fill=Clinical.Stage, y=Proportion, x=Risk.Group)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Clinical Stage for each risk group")+
  scale_fill_manual(values=sample(col_vector, 4))+theme_classic()+xlab("Risk Group")
dev.off()

test_meta$Immune_Subtype = ifelse(test_meta$Immune.Subtype=="C1","Wound Healing",
                             ifelse(test_meta$Immune.Subtype=="C2","IFN-gamma",
                                    ifelse(test_meta$Immune.Subtype=="C3","Inflammatory",
                                           ifelse(test_meta$Immune.Subtype=="C4","Lymphocyte Depleted",
                                                  ifelse(test_meta$Immune.Subtype=="C5","Immunologically Quiet","TGF-beta")))))
df = as.data.frame(t(table(test_meta$Immune_Subtype,test_meta$risk_group)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("TCGA.Subtypes","Risk.Group","Proportion")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_TCGAsubtypes.pdf",
    width=5,height=3)
ggplot(prop, aes(fill=TCGA.Subtypes, y=Proportion, x=Risk.Group)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("TCGA Immune Subtypes for each risk group")+
  scale_fill_manual(values=sample(col_vector, 6))+theme_classic()+xlab("Risk Group")
dev.off()

df = as.data.frame(t(table(test_meta$TCGA.Study,test_meta$risk_group)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("Cancer","Risk.Group","Proportion")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_cancer.pdf",
    width=9,height=6)
ggplot(prop, aes(fill=Cancer, y=Proportion, x=Risk.Group)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Cancer types for each risk group")+
  scale_fill_manual(values=sample(col_vector, 33))+theme_classic()+xlab("Risk Groups")
dev.off()

temp = aggregate(test_meta,by=list(test_meta$risk_group),FUN = mean,na.rm=T)
temp = temp[,c("prolif","exhaustion_full","tissue_memory.full")]
rownames(temp) = temp$Group.1
#temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4")
colnames(temp) = c("Proliferation","Exhaustion","Tissue.Memory")
temp_long = melt(temp)
fun_range <- function(x) {                             
  (x - min(x)) / (max(x) - min(x))
  }
temp_long$Var1 = as.factor(temp_long$Var1)
temp_long$value = fun_range(temp_long$value)
names(temp_long)[1] = "Risk.Groups"
set.seed(66)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_subtype_SigAbundance2.pdf",
    width=7,height=5)
ggplot(data = temp_long,aes(x = Risk.Groups, y = value,fill=Risk.Groups))+
  geom_histogram(stat="identity")+facet_wrap(~Var2,nrow=3)+ 
  theme(text = element_text(size = 25))+scale_fill_manual(values=sample(col_vector, 4))+
  ylab("Score")+xlab("Risk Groups")+theme_classic()
dev.off()

temp = aggregate(test_meta,by=list(test_meta$risk_group),FUN = mean,na.rm=T)
temp = temp[,c("SNV.Neoantigens","Indel.Neoantigens","Silent.Mutation.Rate",
               "Nonsilent.Mutation.Rate","Fraction.Altered","Aneuploidy.Score","Number.of.Segments",
               "Homologous.Recombination.Defects","Intratumor.Heterogeneity")]
rownames(temp) = temp$Group.1
#temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("1","2","3","4")
# temp = temp[,c(1,2,4:10)]
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_subtype_SigAbundance.pdf",
    width=7,height=4)
Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Tumor characteristics across risk groups",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        },row_names_max_width=unit(80,"mm"),column_names_rot = 0)
dev.off()

############### Study figure 6E - TCGA clusters ###############
test_meta_subset = test_meta[-which(test_meta$clinical_stage%in%c("[Not Available]",
                                                                  "[Discrepancy]",
                                                                  "[Not Applicable]","I","IIa","IIb","III","IVa","IVb")),]
test_meta_subset = test_meta[-which(test_meta$Immune.Subtype%in%c("C4","C6")),]
test_meta_subset$Clinical.Stage = ifelse(test_meta_subset$clinical_stage%in%c("Stage I","Stage IA","Stage IA1",
                                                                              "Stage IB","Stage IB1","Stage IB2",
                                                                              "Stage IC"),"Stage I",
                                         ifelse(test_meta_subset$clinical_stage%in%c("Stage II","Stage IIA",
                                                                                     "Stage IIA1","Stage IIA2",
                                                                                     "Stage IIB","Stage IIC"),"Stage II",
                                                ifelse(test_meta_subset$clinical_stage%in%c("Stage III","Stage IIIA",
                                                                                            "Stage IIIB","Stage IIIC",
                                                                                            "Stage IIIC1","Stage IIIC2",
                                                                                            "Stage IS"),"Stage III","Stage IV")))
df = as.data.frame(t(table(test_meta_subset$Clinical.Stage,test_meta_subset$Immune.Subtype)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("Clinical.Stage","TCGA.Subtypes","Proportion")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(6668)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_immunesubtypes_clinicalstage.pdf",
    width=4,height=3)
ggplot(prop, aes(fill=Clinical.Stage, y=Proportion, x=TCGA.Subtypes)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Clinical Stage for each TCGA subgroups")+
  scale_fill_manual(values=sample(col_vector, 4))+theme_classic()+xlab("TCGA Subtypes")
dev.off()

df = as.data.frame(t(table(test_meta_subset$TCGA.Study,test_meta_subset$Immune.Subtype)))
df = reshape2::acast(df, Var2~Var1, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
prop$Var2 = as.factor(prop$Var2)
names(prop) = c("Cancer","TCGA.Subtypes","Proportion")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_immunesubtypes_cancer.pdf",
    ,width=9,height=6)
ggplot(prop, aes(fill=Cancer, y=Proportion, x=TCGA.Subtypes)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Cancer types for each TCGA subtype")+
  scale_fill_manual(values=sample(col_vector, 33))+theme_classic()+xlab("TCGA Subtypes")
dev.off()


test_meta_subset = test_meta[-which(test_meta$Immune.Subtype%in%c("C4","C6")),]
temp = aggregate(test_meta_subset,by=list(test_meta_subset$Immune.Subtype),FUN = mean,na.rm=T)
temp = temp[,c("prolif","exhaustion_full","tissue_memory.full")]
rownames(temp) = temp$Group.1
#temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("C1","C2","C3","C5")
colnames(temp) = c("Proliferation","Exhaustion","Tissue.Memory")
temp_long = melt(temp)
fun_range <- function(x) {                             
  (x - min(x)) / (max(x) - min(x))
}
temp_long$Var1 = as.factor(temp_long$Var1)
temp_long$value = fun_range(temp_long$value)
names(temp_long)[1] = "TCGA.Subtypes"
set.seed(66)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_immunesubtypes_SigAbundance2.pdf",
    width=7,height=5)
ggplot(data = temp_long,aes(x = TCGA.Subtypes, y = value,fill=TCGA.Subtypes))+
  geom_histogram(stat="identity")+facet_wrap(~Var2,nrow=3)+ 
  theme(text = element_text(size = 20))+scale_fill_manual(values=sample(col_vector, 4))+
  ylab("Score")+xlab("TCGA.Subtypes")+theme_classic()
dev.off()

temp = aggregate(test_meta_subset,by=list(test_meta_subset$Immune.Subtype),FUN = mean,na.rm=T)
temp = temp[,c("SNV.Neoantigens","Indel.Neoantigens","Silent.Mutation.Rate",
               "Nonsilent.Mutation.Rate","Fraction.Altered","Aneuploidy.Score","Number.of.Segments",
               "Homologous.Recombination.Defects","Intratumor.Heterogeneity")]
rownames(temp) = temp$Group.1
#temp = temp[,-ncol(temp)]
temp = apply(temp,2,scale)
rownames(temp) = c("C1","C2","C3","C5")
# temp = temp[,c(1,2,4:10)]
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_immunesubtypes_SigAbundance.pdf",
    width=7.5,height=4.5)
Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Tumor characteristics across TCGA subgroups",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        },row_names_max_width=unit(80,"mm"),column_names_rot=0)
dev.off()

ct_fraction = read.delim("/Volumes/she4/hallmark/data/validation-data/TCGA/TCGA.Kallisto.fullIDs.cibersort.relative.tsv")
ct_fraction = ct_fraction[!duplicated(ct_fraction$SampleID),]
ct_fraction$tcga_sample = substr(ct_fraction$SampleID,1,12)
ct_fraction = ct_fraction[!duplicated(ct_fraction$tcga_sample),]
rownames(ct_fraction) = gsub("\\.","-",ct_fraction$tcga_sample)
ct_fraction = ct_fraction[,c(3:24)]
normalize <- function(x, na.rm = TRUE) x*100/sum(x)
ct_fraction = as.data.frame(t(apply(ct_fraction,1,normalize)))
ct_fraction$sampleID = rownames(ct_fraction)
meta_fraction = merge(ct_fraction,train_meta[,c("PI_group","sampleID")],by="sampleID")
meta_fraction = meta_fraction[,-1]
tmp = aggregate(meta_fraction,by=list(meta_fraction$PI_group),FUN=mean)
tmp = tmp[,-c(1,ncol(tmp))]

# meta_fraction = test_meta[,c(205:207,212:239,241)]
# meta_fraction = meta_fraction[complete.cases(meta_fraction),]
# normalize <- function(x, na.rm = TRUE) x*100/sum(x)
# meta_fraction = as.data.frame(t(apply(meta_fraction,1,normalize)))
# tmp = aggregate(meta_fraction,by=list(meta_fraction$PI_group),FUN=mean)
# tmp = tmp[,-c(1,ncol(tmp))]

library(ggplot2)
df = data.frame(group = rep(colnames(tmp),4),
                cluster = rep(c(1:4),each=22),
                value = c(as.vector(unlist(tmp[1,])),
                          as.vector(unlist(tmp[2,])),
                          as.vector(unlist(tmp[3,])),
                          as.vector(unlist(tmp[4,]))))
n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
df$cluster = as.factor(df$cluster)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_subtype_celltypecomposition.pdf",
    width=8,height=5)
ggplot(df, aes(x=cluster, y=value, fill=group))+
  geom_bar(stat = "identity")+facet_wrap(~group)+
  scale_fill_manual(values = col_vector)
dev.off()

ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(stat = "identity")+facet_wrap(~cluster)+
  scale_fill_manual(values = col_vector)+coord_polar("y", start=0)

# temp = aggregate(test_meta,by=list(train_meta$PI_group),FUN = mean,na.rm=T)
# temp = temp[,213:240]
# fun_range <- function(x) {                              # Create user-defined function
#   (x - min(x)) / (max(x) - min(x))
# }
# normalize <- function(x, na.rm = TRUE) x*100/sum(x)
# temp = as.data.frame(t(apply(temp,1,fun_range)))
# temp = as.data.frame(t(apply(temp,1,normalize)))
# df = data.frame(group = rep(colnames(temp),4),
#                 cluster = rep(c(1:4),each=28),
#                 value = c(as.vector(unlist(temp[1,])),
#                           as.vector(unlist(temp[2,])),
#                           as.vector(unlist(temp[3,])),
#                           as.vector(unlist(temp[4,]))))
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/TCGA_PI_subtype_celltypecomposition.pdf",
#     width=10,height=8)
# ggplot(df, aes(x="", y=value, fill=group))+
#   geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+facet_wrap(~cluster)
# dev.off()

pro = c("IFNG","IL1A","IL1B","IL2")
anti = c("IL4","IL10","IL11","TGFB1")
pro_anti = list(pro,anti)
names(pro_anti) = c("Pro-inflammatory","Anti-inflammatory")
pdf("~/desktop/umap_kmean.pdf",width=6,height=5)
ggplot(dat,aes(x=V1,y=V2,colour=type))+geom_point()+xlab("UMAP 1")+ylab("UMAP 2")
dev.off()

test <- data.frame(Truth = full_data$Response_binary, 
                   MP=pred_MP,Lymphoid_MP=pred_L,
                   Myeloid_MP=pred_M,
                   Hallmark=pred_H,
                   KEGG=pred_K,stringsAsFactors = FALSE)
test = melt(test, id.vars="Truth")
names(test)[2] = "Gene_sets"
test$Gene_sets = ifelse(test$Gene_sets=="MP","MP (AUC: 0.85)",
                        ifelse(test$Gene_sets=="Lymphoid_MP","Lymphoid_MP (AUC: 0.73)",
                               ifelse(test$Gene_sets=="Myeloid_MP","Myeloid_MP (AUC: 0.70)",
                                      ifelse(test$Gene_sets=="Hallmark","Hallmark (AUC: 0.87)","KEGG (AUC: 0.72)"))))
test$Gene_sets = factor(test$Gene_sets,levels=c("Hallmark (AUC: 0.87)","MP (AUC: 0.85)",
                                                "Lymphoid_MP (AUC: 0.73)","Myeloid_MP (AUC: 0.70)",
                                                "KEGG (AUC: 0.72)"))
names(test)[2] = "Gene.Sets"
pdf("~/desktop/ICB_ROC.pdf",width=6,height=4)
ggplot(test, aes(d = Truth, color = Gene.Sets, m=value)) + 
  geom_roc(n.cuts = 0)+ggtitle("ICB Response classification accuracy")+
  theme(text = element_text(size = 15))+theme_classic()+
  annotate(geom="text", x=0.4, y=0.5, 
           label="MP vs Hallmark: p.value=0.31",color="black",size = 3,hjust=0)+
  annotate(geom="text", x=0.4, y=0.45, 
           label="MP vs KEGG: p.value=0.07",color="black",size = 3,hjust=0)
dev.off()
# roc1 <- roc(full_data$Response_binary,pred_H)
# roc2 <- roc(full_data$Response_binary,pred_MP)
# roc3 <- roc(full_data$Response_binary,pred_M)
# roc4 <- roc(full_data$Response_binary,pred_L)
# roc5 <- roc(full_data$Response_binary,pred_K)
# roc.test(roc2, roc1,alternative="two.sided",method="bootstrap")
# roc.test(roc2, roc5,alternative="two.sided",method="bootstrap")


test <- data.frame(Truth = full_data$Response_binary, 
                   MP=pred_mp,Variables=pred_v,Both=pred,
                   stringsAsFactors = FALSE)
test = melt(test, id.vars="Truth")
test$Gene_sets = ifelse(test$Gene_sets=="MP","MP (AUC: 0.85)",
                        ifelse(test$Gene_sets=="Variables","Variables (AUC: 0.58)", "Both (AUC: 0.92)"))
test$Gene_sets = factor(test$Gene_sets,levels=c("Both (AUC: 0.92)","MP (AUC: 0.85)","Variables (AUC: 0.58)"))
names(test)[2] = "Gene.Sets"
pdf("~/desktop/ICB_ROC2.pdf",width=6,height=4)
ggplot(test, aes(d = Truth, color = Gene.Sets, m=value)) + 
  geom_roc(n.cuts = 0)+ggtitle("ICB Response classification accuracy")+
  theme(text = element_text(size = 15))+theme_classic()+
  annotate(geom="text", x=0.4, y=0.5, 
           label="MP vs Variables: p.value<0.01",color="black",size = 3,hjust=0)+
  annotate(geom="text", x=0.4, y=0.45, 
           label="Both vs Variables: p.value<0.0001",color="black",size = 3,hjust=0)
dev.off()
roc1 <- roc(full_data$Response_binary,pred_mp)
roc2 <- roc(full_data$Response_binary,pred_v)
roc3 <- roc(full_data$Response_binary,pred)
roc.test(roc1, roc2,alternative="two.sided",method="bootstrap")
roc.test(roc3, roc2,alternative="two.sided",method="bootstrap")

fit<- survfit(Surv(meta$OS.time, meta$OS.x) ~ Immune.Subtype, data = meta)
pdf("~/desktop/TCGA_KM_plot.pdf",width=7,height=5,onefile = F)
ggsurvplot(fit, data = meta,pval=TRUE,pval.method=T)+ggtitle("Kaplan Meier for TCGA Pan-cancer")
dev.off()

survdiff(Surv(meta$OS.time, meta$OS.x) ~ Immune.Subtype, data = meta)

Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Signature expression across different immune subtypes",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        })

table(meta$Immune.Subtype_new)
count[1:5,1:5]

gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
H_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(H_list) = unique(gene_sets$gs_name)

kegg_sets = read.csv(paste0(path_data,"kegg_immune_sets.csv"))[,-1]
KEGG_list = list()
for (i in 1:length(unique(kegg_sets$name))){
  KEGG_list[[i]] = unique(kegg_sets[which(kegg_sets$name==unique(kegg_sets$name)[i]),"symbol"])
}
names(KEGG_list) = unique(kegg_sets$name)

H = corto::ssgsea(inmat = as.matrix(count),groups = H_list)
H = as.data.frame(t(H))
KEGG = corto::ssgsea(inmat = as.matrix(count),groups = KEGG_list)
KEGG = as.data.frame(t(KEGG))

KEGG$Immune.Subtype_new = as.factor(kmean_cluster$cluster)
temp = aggregate(KEGG,by=list(KEGG$Immune.Subtype_new),FUN = mean)
temp = temp[,c(1:30)]
rownames(temp) = temp$Group.1
temp = temp[,-1]
names(temp) = names(KEGG_list)[1:29]
pdf("~/desktop/TCGA_novel_subtype_abundance.pdf",width=12,height=8)
Heatmap(as.matrix(t(temp)),
        name = "Expression",
        cluster_rows = T,cluster_columns = F,show_column_dend = F,show_row_dend = F,
        column_title = "Signature expression across different immune subtypes",
        column_names_max_height=unit(150,"mm"),cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", as.matrix(t(temp))[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


