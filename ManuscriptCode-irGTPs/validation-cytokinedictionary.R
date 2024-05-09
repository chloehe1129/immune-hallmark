library(dplyr)
library(Seurat)
library(fgsea)
library(ComplexHeatmap)
library(clusterProfiler)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

annotation_L = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
annotation_M = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))
annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
gs_list_L = list()
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names)){
  gs_list_L[[i]] = read.csv(file_names[i],row.names = 1)[,1]
}
names(gs_list_L) = names(annotation)
annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))
gs_list_M = list()
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names)){
  gs_list_M[[i]] = read.csv(file_names[i],row.names = 1)[,1]
}
names(gs_list_M) = names(annotation)
names(annotation_M)[4] = paste0(names(annotation_M)[4],"_M")

TERM2GENE = NULL
#file_names_MP <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = "_based")
file_names_MP <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv"))
#where you have your files
for (i in 1:length(file_names_MP)){
  TERM2GENE = rbind(TERM2GENE,cbind(rep(names(annotation_M)[i],50),read.csv(file_names_MP[i])[,2]))
}
TERM2GENE = as.data.frame(TERM2GENE)
names(TERM2GENE) = c("term",'gene')
TERM2GENE_L = NULL
file_names_MP <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names_MP)){
  TERM2GENE_L = rbind(TERM2GENE_L,cbind(rep(names(annotation_L)[i],50),read.csv(file_names_MP[i])[,2]))
}
TERM2GENE_L = as.data.frame(TERM2GENE_L)
names(TERM2GENE_L) = c("term",'gene')
TERM2GENE = rbind(TERM2GENE,TERM2GENE_L)

all_dic = list.files(paste0(path_data,"cytokine-dictionary/"),full.names = T,pattern="human")
dic = readRDS(all_dic[1])
for(i in 2:length(all_dic)){
  dic = merge(dic,readRDS(all_dic[i]))
  print(i)
}
dic = NormalizeData(dic)
saveRDS(dic,paste0(path_data,"cytokine-dictionary/all_dic.rds"))

dic = readRDS(paste0(path_data,"cytokine-dictionary/all_dic.rds"))
cytokines = setdiff(unique(dic$sample),"PBS")

mp = 23

NES = matrix(0,86,15)
pvalue = matrix(0,86,15)
rownames(NES) = rownames(pvalue) = cytokines
colnames(NES) = colnames(pvalue) = unique(dic$celltype)

Idents(dic) = "celltype"
for(ct in 1:length(unique(dic$celltype))){
  for(cyto in cytokines){
    if(length(which(subset(dic,subset=celltype==unique(dic$celltype)[ct])$sample==cyto))<3){
      NES[cyto,ct] = 0
      pvalue[cyto,ct] = 1
    }else{
      markers <- FindMarkers(dic, ident.1 = cyto, ident.2 = "PBS",group.by = "sample",
                             subset.ident = unique(dic$celltype)[ct],
                             logfc.threshold = 0,min.diff.pct=0,min.pct=0)
      markers = markers[order(markers$avg_log2FC,decreasing = T),]
      gene_list = markers$avg_log2FC
      names(gene_list) = rownames(markers)
      gsea = GSEA(gene_list,TERM2GENE = TERM2GENE[which(TERM2GENE$term==unique(TERM2GENE$term)[mp]),],pvalueCutoff=1,minGSSize=0)
      if(dim(gsea@result)[1]==0){
        NES[cyto,ct] = 0
        pvalue[cyto,ct] = 1
      }else{
        NES[cyto,ct] = gsea@result$NES
        pvalue[cyto,ct] = gsea@result$p.adjust
      }
    }
    print(cyto)
  }
  print(ct)
  write.csv(NES,paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],".csv"))
  write.csv(pvalue,paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],"_pvalue.csv"))
}
# if(length(which(rownames(NES)%in%na_cyto))!=0){
#   NES = NES[-which(rownames(NES)%in%na_cyto),]
# }
mp=23
wanted = c("B_cell","T_cell_CD4","T_cell_CD8","Treg","Macrophage","cDC1","NK_cell")
NES = read.csv(paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],".csv"),row.names = 1)
pvalue = read.csv(paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],"_pvalue.csv"),row.names = 1)
NES = NES[,wanted]
pvalue = pvalue[,wanted]
# pdf(paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],".pdf"),height=16,width = 12)
# print(Heatmap(NES,name="NES",column_title = unique(TERM2GENE$term)[mp],na_col = "grey",
#               cluster_rows = T,cluster_columns = T,
#               column_names_max_height = unit(9, "cm"),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#               }))
# dev.off()
df_IFN = as.data.frame(NES[grep("IFN",rownames(NES)),])
df_IFN$Cytokine = "Type I Interferon"
df_other = as.data.frame(NES[-grep("IFN",rownames(NES)),])
df_other$Cytokine = "Non-Interferon"
df = as.data.frame(rbind(df_IFN,df_other))
# Heatmap(t(df[,1:7]),cluster_columns = F)
df_long = melt(setDT(df),id.vars = "Cytokine")
df_long$variable = factor(df_long$variable,
                          levels = c("cDC1","cDC2","eTAC","ILC","Macrophage","MigDC",
                                     "Monocyte","Neutrophil","NK_cell","pDC",
                                     "B_cell","T_cell_CD4","T_cell_CD8","T_cell_gd","Treg"))
aggregate(df_long,by=list(df_long$Cytokine,df_long$variable),FUN=mean)
ggplot(data=df_long,aes(x=value,fill=Cytokine))+geom_boxplot(outlier.shape = NA)+
  facet_wrap(~variable,scale='free_y',ncol=1)+ylab(unique(TERM2GENE$term)[mp])+xlab("Gene set activity")+
  theme_classic()+
  theme(legend.position = 'bottom', legend.direction = "horizontal")
pdf(paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],".pdf"),
    height=6,width = 5)
ggplot(data=df_long,aes(x=value,fill=Cytokine))+geom_density(alpha=0.4)+
  facet_wrap(~variable,scale='free_y',ncol=1)+ylab(unique(TERM2GENE$term)[mp])+xlab("Gene set activity")+
  theme_classic()+
  theme(legend.position = 'bottom', legend.direction = "horizontal")
  # stat_compare_means(method="t.test",label="p.signif")+
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
dev.off()
# for(mp in 1:28){
#   
#   NES = matrix(0,86,15)
#   pvalue = matrix(0,86,15)
#   rownames(NES) = rownames(pvalue) = cytokines
#   colnames(NES) = colnames(pvalue) = unique(dic$celltype)
#     
#   Idents(dic) = "celltype"
#   for(ct in 1:length(unique(dic$celltype))){
#     for(cyto in cytokines){
#       if(length(which(subset(dic,subset=celltype==unique(dic$celltype)[ct])$sample==cyto))<3){
#         NES[cyto,ct] = NA
#         pvalue[cyto,ct] = NA
#         }else{
#           markers <- FindMarkers(dic, ident.1 = cyto, ident.2 = "PBS",group.by = "sample",subset.ident = unique(dic$celltype)[ct])
#           markers = markers[order(markers$avg_log2FC,decreasing = T),]
#           gene_list = markers$avg_log2FC
#           names(gene_list) = rownames(markers)
#           gsea = GSEA(gene_list,TERM2GENE = TERM2GENE[which(TERM2GENE$term==unique(TERM2GENE$term)[2]),],pvalueCutoff=1,minGSSize=0)
#           if(dim(gsea@result)[1]==0){
#             NES[cyto,ct] = NA
#             pvalue[cyto,ct] = NA
#           }else{
#             NES[cyto,ct] = gsea@result$NES
#             pvalue[cyto,ct] = gsea@result$p.adjust
#           }
#         }
#       print(cyto)
#     }
#     print(ct)
#   }
#   # if(length(which(rownames(NES)%in%na_cyto))!=0){
#   #   NES = NES[-which(rownames(NES)%in%na_cyto),]
#   # }
#   pdf(paste0(path_result,"validation/cytokine-dictionary-",unique(TERM2GENE$term)[mp],".pdf"),height=16,width = 12)
#   print(Heatmap(NES,name="NES",column_title = unique(TERM2GENE$term)[mp],na_col = "grey",
#                 cluster_rows = F,cluster_columns = F,
#                 column_names_max_height = unit(9, "cm")))
#   dev.off()
#   print(mp)
# }

# for(i in 1:length(all_dic)){
#   
#   NES = matrix(0,86,28)
#   pvalue = matrix(0,86,28)
#   rownames(NES) = rownames(pvalue) = cytokines
#   
#   dic = readRDS(all_dic[i])
#   dic <- NormalizeData(dic)
#   Idents(dic) = "sample"
#   cytokines = setdiff(unique(dic$sample),"PBS")
#   
#   na_cyto = NULL
#   for(cyto in cytokines){
#     if(length(which(dic$sample==cyto))<3){
#       NES[cyto,] = NA
#       pvalue[cyto,] = NA
#       na_cyto = c(na_cyto,cyto)
#     }else{
#       markers <- FindMarkers(dic, ident.1 = cyto, ident.2 = "PBS")
#       markers = markers[order(markers$avg_log2FC,decreasing = T),]
#       gene_list = markers$avg_log2FC
#       names(gene_list) = rownames(markers)
#       gsea = GSEA(gene_list,TERM2GENE = TERM2GENE,pvalueCutoff=1,minGSSize=0)
#       
#       NES[cyto,] = gsea@result[unique(TERM2GENE$term),"NES"]
#       pvalue[cyto,] = gsea@result[unique(TERM2GENE$term),"p.adjust"]
#     }
#     
#     # gsea = gsea[order(gsea$NES,decreasing = T),]
#     # gsea$Description <- factor(gsea$Description, levels = rev(rownames(gsea)))
#     # 
#     # print(ggplot(gsea, aes(x = NES, y = Description)) +
#     #   geom_point(aes(size = NES, color = p.adjust)) +
#     #   theme_bw(base_size = 14) +
#     #   scale_colour_gradient(limits=c(0, 1), low="brown2",high="cornflowerblue") +
#     #   ylab(NULL) +xlab("Normalized Enrichment Score")+
#     #   ggtitle(paste0("Enrichment for ",cyto,"-simulated cells"))+
#     #   theme(text=element_text(size=13,face = "bold")))
#   }
#   
#   colnames(NES) = unique(TERM2GENE$term)
#   if(length(which(rownames(NES)%in%na_cyto))!=0){
#     NES = NES[-which(rownames(NES)%in%na_cyto),]
#   }
#   pdf(paste0(path_result,"validation/cytokine-dictionary-",substr(all_dic[i],65,nchar(all_dic[i])-4),".pdf"),height=16,width = 12)
#   print(Heatmap(NES,name="NES",column_title = substr(all_dic[i],65,nchar(all_dic[i])-4),na_col = "grey",
#           column_names_max_height = unit(9, "cm")))
#   dev.off()
# 
# }

