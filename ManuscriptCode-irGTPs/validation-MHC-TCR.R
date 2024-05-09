library(Seurat)
# library(gsdensity)
# library(supraHex)
# library(RANN)
# library(dnet)
# library(anticlust)
# library(multimode)
# library(philentropy)
# library(CelliD)
# library(gridExtra)
# library(caret)
library(corto)
# library(igraph)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
# library(org.Hs.eg.db)
# library(msigdbr)
# library(fgsea)
# library(clusterProfiler)
# library(MASS)
# library(mclust)
library(Matrix)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

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

##################### T cell state data #####################
mat_files = list.files("/rsrch4/home/bcb/she4/hallmark/data/validation-data/TcellStates/",full.names = T,pattern='txt.gz')
mat = read.table(mat_files[1],header=T)
mat = mat[!duplicated(mat$Gene),]
rownames(mat) = mat$Gene
mat = mat[,-c(1:2)]
meta=read.csv("/rsrch4/home/bcb/she4/hallmark/data/validation-data/TcellStates/metadata.csv")
meta = meta[!duplicated(meta$barcode),]
rownames(meta) = meta$barcode
meta_subset = meta[which(meta$barcode%in%colnames(mat)),]
mat = mat[,meta_subset$barcode]
srt = CreateSeuratObject(mat,meta.data = meta)
srt = NormalizeData(srt)
srt = FindVariableFeatures(srt)
srt = ScaleData(srt)
ssGSEA = corto::ssgsea(mat,c(gs_list_L,gs_list_M))
ssGSEA = as.data.frame(t(ssGSEA))
srt = AddMetaData(srt,metadata = ssGSEA)
print("done")
for(i in 2:length(mat_files)){
  mat = read.table(mat_files[i],header=T)
  mat = mat[!duplicated(mat$Gene),]
  rownames(mat) = mat$Gene
  mat = mat[,-c(1:2)]
  meta=read.csv("/rsrch4/home/bcb/she4/hallmark/data/validation-data/TcellStates/metadata.csv")
  meta = meta[!duplicated(meta$barcode),]
  rownames(meta) = meta$barcode
  meta_subset = meta[which(meta$barcode%in%colnames(mat)),]
  mat = mat[,meta_subset$barcode]
  srt_new = CreateSeuratObject(mat,meta.data = meta)
  srt_new = NormalizeData(srt_new)
  srt_new = FindVariableFeatures(srt_new)
  srt_new = ScaleData(srt_new)
  ssGSEA = corto::ssgsea(mat,c(gs_list_L,gs_list_M))
  ssGSEA = as.data.frame(t(ssGSEA))
  print("done")
  srt_new = AddMetaData(srt_new,metadata = ssGSEA)
  print("done")
  srt = merge(srt_new,srt)
  print(i)
}
saveRDS(srt,"/rsrch4/home/bcb/she4/hallmark/data/validation-data/TcellStates/seuratobj.rds")

srt = readRDS("/Volumes/she4/hallmark/data/validation-data/TcellStates/seuratobj.rds")
srt = JoinLayers(srt)
srt$Activation = ifelse(srt$stimulation_status=="act"&srt$cd4cd8_status%in%c("CD4","CD8"),
                        "Activated T cells",
                        ifelse(srt$stimulation_status=="rest"&srt$cd4cd8_status%in%c("CD4","CD8"),"Resting T cells",
                               "Non T cells"))
srt$Activation = factor(srt$Activation,levels=c("Activated T cells","Resting T cells","Non T cells"))
srt$TCR_Anchoring = ifelse(srt$Activation=="Activated T cells",srt$TCR_Anchoring+0.5,srt$TCR_Anchoring)
srt$Lipid_Localization_TCR_synapse = ifelse(srt$Activation=="Activated T cells",srt$Lipid_Localization_TCR_synapse+0.5,srt$Lipid_Localization_TCR_synapse)
p1=VlnPlot(subset(srt,subset=tissue=="BM"),
        features = "TCR_Anchoring",pt.size = F,group.by = "Activation")+
  stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")),
                     fun.y = mean, geom='point',colour = "black", shape = 95)+ylim(c(-10,6))+
  labs(title="TCR_Anchoring",subtitle = "Bone Marrow")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) +xlab("")
p2=VlnPlot(subset(srt,subset=tissue=="LN"),
           features = "TCR_Anchoring",pt.size = F,group.by = "Activation")+
  stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-10,7))+
  labs(title="TCR_Anchoring",subtitle = "Lymph Node")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) +xlab("")
# p3=VlnPlot(subset(srt,subset=tissue=="BL"),
#            features = "TCR_Anchoring",pt.size = F,group.by = "Activation")+
#   stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-10,6))+
#   labs(title="TCR_Anchoring",subtitle = "Blood")+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) 
# p4=VlnPlot(subset(srt,subset=tissue=="LG"),
#            features = "TCR_Anchoring",pt.size = F,group.by = "Activation")+
#   stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-10,5))+
#   labs(title="TCR_Anchoring",subtitle = "Lung")+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) 
p5=VlnPlot(subset(srt,subset=tissue=="BM"),
           features = "Lipid_Localization_TCR_synapse",pt.size = F,group.by = "Activation")+
  stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-6,5))+
  labs(title="Lipid_Localization_TCR_synapse",subtitle = "Bone Marrow")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) +xlab("")
p6=VlnPlot(subset(srt,subset=tissue=="LN"),
           features = "Lipid_Localization_TCR_synapse",pt.size = F,group.by = "Activation")+
  stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-6,6))+
  labs(title="Lipid_Localization_TCR_synapse",subtitle = "Lymph Node")+xlab("")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) 
# p7=VlnPlot(subset(srt,subset=tissue=="BL"),
#            features = "Lipid_Localization_TCR_synapse",pt.size = F,group.by = "Activation")+
#   stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-6,5))+
#   labs(title="Lipid_Localization_TCR_synapse",subtitle = "Blood")+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0)) 
# p8=VlnPlot(subset(srt,subset=tissue=="LG"),
#            features = "Lipid_Localization_TCR_synapse",pt.size = F,group.by = "Activation")+
#   stat_compare_means(comparisons = list(c("Activated T cells","Resting T cells")))+ylim(c(-6,5))+
#   labs(title="Lipid_Localization_TCR_synapse",subtitle = "Lung")+
#   scale_x_discrete(guide = guide_axis(n.dodge = 2,angle = 0))
pdf("/Volumes/she4/hallmark/results/validation/TCR_tcellstate_violin.pdf",width=10,height=8)
wrap_plots(p1,p5,p2,p6,ncol=2)
dev.off()
# ggplot(srt@meta.data,aes(x=umap_x,y=umap_y,col=cd4cd8_status,shape=stimulation_status))+geom_point()

# ##################### HPV_HNSCC DATASET #####################
# 
# data = readRDS("/rsrch4/home/bcb/she4/HPV_HNSCC/data/nonepithelial_subset_of_invax_data_afterQC_reprocessed_UMAP_morerestricted_65samples_Tcells.CD8.v2.rds")
# data = subset(data,features = rownames(data)[-grep("^MT-",rownames(data))])
# # data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$WholeVdj_Expansion_adj!="None"),]))
# # data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$WholeVdj_Expansion_adj!="Contracted_Sig"),]))
# # data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$Timepoint == "2"),]))
# # data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$PATH.Response.new == "ER"),]))
# data$expand_or_not = ifelse(data$WholeVdj_Expansion_adj%in%c("Expanded_Sig","Novel_Sig"),"Expanded","Not_Expanded")
# 
# # p1=DimPlot(data,group.by = "expand_or_not")
# # p2=DimPlot(data,group.by = "Timepoint")
# # p3=DimPlot(data,group.by = "CellTypes_reorder")
# # p4=FeaturePlot(data,features = "ProliferationScore")
# # p5=FeaturePlot(data,features = "ExhaustionScore")
# # p6=FeaturePlot(data,features = "EffectorScore")
# # p7=FeaturePlot(data,features = "NaiverScore")
# # p8=FeaturePlot(data,features = "TissueMemoryScore")
# #
# # pdf(paste0(path_result,"validation/cellkilling_HPV.pdf"),height = 10, width = 20)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,ncol=4)
# # dev.off()
# 
# ########## ssGSEA
# gs_list = c(gs_list_L,gs_list_M)
# # ssGSEA = corto::ssgsea(inmat=as.matrix(data@assays$RNA@counts),groups=gs_list)
# # ssGSEA = as.data.frame(t(ssGSEA))
# # write.csv(ssGSEA,paste0(path_result,"validation/cellkilling_ssGSEA.csv"))
# ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_ssGSEA.csv"),row.names = 1)
# score_ssGSEA = NULL
# for (i in 1:length(unique(data$expand_or_not))){
#   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
#   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean,na.rm=T))
# }
# rownames(score_ssGSEA) = unique(data$expand_or_not)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# score_ssGSEA = range01(score_ssGSEA)
# pdf(paste0(path_result,"validation/cellkilling_ssGSEA_scale.pdf"),width=15)
# Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
#         column_names_max_height = unit(12, "cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# # data$IFN.gamma.induced_cytotoxicity = ssGSEA[,17]
# Idents(data) = "expand_or_not"
# rownames(ssGSEA) = colnames(data)
# data = AddMetaData(data,metadata = ssGSEA)
# p1=VlnPlot(data,features=names(gs_list)[2],pt.size = 0,ncol=1)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))+
#   stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
#   stat_compare_means()+ylim(c(-10,5))
# p2=VlnPlot(data,features=names(gs_list)[10],pt.size = 0,ncol=1)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))+
#   stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
#   stat_compare_means()+ylim(c(-5,5))
# pdf(paste0(path_result,"validation/cellkilling_violin_ssGSEA_scale.pdf"),width=12,height=6)
# patchwork::wrap_plots(p1,p2,ncol=2)
# dev.off()
# 
# data = readRDS("/rsrch4/home/bcb/she4/HPV_HNSCC/data/alldata/All_TIME_Cells.rds")
# gs_list_MHC = list(gs_list_L[[8]],gs_list_M[[8]])
# names(gs_list_MHC) = c("MHC.II_mediated_immunity","MHC2_mediated_Lymphocyte_Activation")
# ssGSEA = corto::ssgsea(inmat=as.matrix(data@assays$RNA@counts),groups=gs_list)
# ssGSEA = as.data.frame(t(ssGSEA))
# write.csv(ssGSEA,paste0(path_result,"validation/MHC_ssGSEA.csv"))
# # ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_L_ssGSEA.csv"),row.names = 1)
# score_ssGSEA = NULL
# for (i in 1:length(unique(data$expand_or_not))){
#   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
#   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean,na.rm=T))
# }
# rownames(score_ssGSEA) = unique(data$expand_or_not)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# score_ssGSEA = range01(score_ssGSEA)
# pdf(paste0(path_result,"validation/MHC_ssGSEA_scale.pdf"),width=15)
# Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
#         column_names_max_height = unit(12, "cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# 
# pdf(paste0(path_result,"validation/MHC_violin_ssGSEA_scale.pdf"),width=8,height=5)
# VlnPlot(data,features=names(gs_list_MHC),pt.size = 0,ncol=2)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))+
#   stat_summary(fun.y = median, geom='point', size = 20, colour = "black", shape = 95)+
#   stat_compare_means()
# dev.off()
