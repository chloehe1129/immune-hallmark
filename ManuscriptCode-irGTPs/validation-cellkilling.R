library(Seurat)
library(gsdensity)
library(supraHex)
library(RANN)
library(dnet)
library(anticlust)
library(multimode)
library(philentropy)
library(CelliD)
library(gridExtra)
library(caret)
library(corto)
library(igraph)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(MASS)
library(mclust)

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

##################### HPV_HNSCC DATASET #####################

# data = readRDS("/rsrch4/home/bcb/she4/HPV_HNSCC/data/nonepithelial_subset_of_invax_data_afterQC_reprocessed_UMAP_morerestricted_65samples_Tcells.CD8.v2.rds")
# data = subset(data,features = rownames(data)[-grep("^MT-",rownames(data))])
# data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$WholeVdj_Expansion_adj!="None"),]))
# data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$WholeVdj_Expansion_adj!="Contracted_Sig"),]))
# data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$Timepoint == "2"),]))
# data = subset(data,cells=rownames(data@meta.data[which(data@meta.data$PATH.Response.new == "ER"),]))
# data$expand_or_not = ifelse(data$WholeVdj_Expansion_adj%in%"Expanded_Sig","Expanded","Not_Expanded")
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
########## ssGSEA

# ssGSEA = corto::ssgsea(inmat=as.matrix(data@assays$RNA@counts),groups=gs_list_L)
# ssGSEA = as.data.frame(t(ssGSEA))
# write.csv(ssGSEA,paste0(path_result,"validation/cellkilling_L_ssGSEA.csv"))
# # ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_L_ssGSEA.csv"),row.names = 1)
# score_ssGSEA = NULL
# for (i in 1:length(unique(data$expand_or_not))){
#   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
#   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean,na.rm=T))
# }
# rownames(score_ssGSEA) = unique(data$expand_or_not)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# score_ssGSEA = range01(score_ssGSEA)
# pdf(paste0(path_result,"validation/cellkilling_L_ssGSEA_scale.pdf"),width=15)
# Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
#         column_names_max_height = unit(12, "cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# data$IFN.gamma.induced_cytotoxicity = ssGSEA[,17]
# Idents(data) = "expand_or_not"
# pdf(paste0(path_result,"validation/cellkilling_violin_ssGSEA_scale.pdf"),width=6,height=6)
# VlnPlot(data,features="IFN.gamma.induced_cytotoxicity",pt.size = 0)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))+
#   stat_summary(fun.y = median, geom='point', size = 20, colour = "black", shape = 95)+
#   stat_compare_means()
# dev.off()

# ssGSEA = corto::ssgsea(inmat=as.matrix(data@assays$RNA@counts),groups=gs_list_M)
# ssGSEA = as.data.frame(t(ssGSEA))
# # write.csv(ssGSEA,paste0(path_result,"validation/cellkilling_M_ssGSEA.csv"))
# # ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_M_ssGSEA.csv"),row.names = 1)
# score_ssGSEA = NULL
# for (i in 1:length(unique(data$expand_or_not))){
#   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
#   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean))
# }
# score_ssGSEA = score_ssGSEA
# rownames(score_ssGSEA) = unique(data$expand_or_not)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# score_ssGSEA = range01(score_ssGSEA)
# pdf(paste0(path_result,"validation/cellkilling_M_ssGSEA_scale.pdf"),width=13)
# Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
#         column_names_max_height = unit(12, "cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# #
# # ########## GSD
# #
# # ce <- compute.mca(object = data)
# # res <- compute.kld(coembed = ce,
# #                    genes.use = intersect(rownames(ce), rownames(data)), # this intersection is to select only genes, not cells.
# #                    n.grids = 100,
# #                    gene.set.list = c(gs_list_L,gs_list_M),
# #                    gene.set.cutoff = 3,
# #                    n.times = 100)
# # cells <- colnames(data)
# #
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# #
# # res = matrix(0,ncol(data),19)
# # for (gs in 1:19){
# #   res[,gs] <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
# # }
# # write.csv(res,paste0(path_result,"validation/cellkilling_L_GSD_scale.csv"))
# #
# # res = matrix(0,ncol(data),9)
# # for (gs in 1:9){
# #   res[,gs] <- run.rwr(el = el, gene_set = gs_list_M[[gs]], cells = cells)
# # }
# # write.csv(res,paste0(path_result,"validation/cellkilling_M_GSD_scale.csv"))
# #
# # ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_L_GSD_scale.csv"),row.names = 1)
# # score_ssGSEA = NULL
# # for (i in 1:length(unique(data$expand_or_not))){
# #   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
# #   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean))
# # }
# # #score_ssGSEA = score_ssGSEA[-2,]
# # rownames(score_ssGSEA) = unique(data$expand_or_not)
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # score_ssGSEA = range01(score_ssGSEA)
# # pdf(paste0(path_result,"validation/cellkilling_L_GSD_scale.pdf"),width=15)
# # Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
# #         column_names_max_height = unit(12, "cm"),
# #         cell_fun = function(j, i, x, y, width, height, fill) {
# #           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
# #         })
# # dev.off()
# #
# # ssGSEA = read.csv(paste0(path_result,"validation/cellkilling_M_GSD_scale.csv"),row.names = 1)
# # score_ssGSEA = NULL
# # for (i in 1:length(unique(data$expand_or_not))){
# #   score_subset = ssGSEA[which(data$expand_or_not==unique(data$expand_or_not)[i]),]
# #   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean))
# # }
# # score_ssGSEA = score_ssGSEA[-2,]
# # rownames(score_ssGSEA) = unique(data$expand_or_not)[-2]
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # score_ssGSEA = range01(score_ssGSEA)
# # pdf(paste0(path_result,"validation/cellkilling_M_GSD_scale.pdf"),width=15)
# # Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "HPV_HNSCC",
# #         column_names_max_height = unit(12, "cm"),
# #         cell_fun = function(j, i, x, y, width, height, fill) {
# #           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
# #         })
# # dev.off()
# 
# ##################### Cell Killing/Exhaustion/Others with COVID-19 data #####################
# 
srt <- readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu.rds"))
srt = subset(srt,cells = colnames(srt)[which(srt@meta.data$cell.type%in%c("B","DC","CD4 T","CD4m T","CD4n T","CD8eff T","CD8m T","gd T","NK"))])
srt@meta.data$drug = ifelse(srt@meta.data$Donor%in%c("C1","C3","C4"),"Azithromycin",ifelse(srt@meta.data$Donor%in%c("H1","H2","H3","H4","H5","H6"),"Healthy","Remdesivir"))
Idents(srt) = "cell.type.fine"

ce <- compute.mca(object = srt)
cells <- colnames(srt)
source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
el <- compute.nn.edges(coembed = ce, nn.use = 3000)

for (gs in 1:20){
  if(gs != 20){
    cv <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
    srt@meta.data$geneset <- cv[colnames(srt)]
    colnames(srt@meta.data)[which(colnames(srt@meta.data)=="geneset")] = names(gs_list_L[gs])
    assign(paste0("p",gs),FeaturePlot(srt,features = names(gs_list_L[gs]),
                                      label=T,raster = T)+ggtitle(names(gs_list_L[gs]))+ NoLegend()+
             theme(plot.title = element_text(face = "bold",size=15)))
  }else{
    cv <- run.rwr(el = el, gene_set = c("PRF1","GZMA","GZMB","GZMH","GNLY"), cells = cells)
    srt@meta.data$geneset <- cv[colnames(srt)]
    colnames(srt@meta.data)[which(colnames(srt@meta.data)=="geneset")] = "Cytotoxicity.Score"
  }
}
saveRDS(srt,paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD.rds"))
srt = readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD.rds"))
for (gs in c(1,6,11)){
  cv <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
  srt@meta.data$geneset <- cv[colnames(srt)]
  colnames(srt@meta.data)[which(colnames(srt@meta.data)=="geneset")] = names(gs_list_L[gs])
  # assign(paste0("p",gs),FeaturePlot(srt,features = names(gs_list_L[gs]),
  #                                   label=T,raster = T)+ggtitle(names(gs_list_L[gs]))+ NoLegend()+
  #          theme(plot.title = element_text(face = "bold",size=15)))
}
saveRDS(srt,paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD-L.rds"))

srt = readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD-L.rds"))
pdf(paste0(path_result,"validation/CellCycle_Lymphocyte_all_COVID19_UMAP.pdf"),height = 26, width = 32)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,nrow=4,ncol=5)
dev.off()

p1=VlnPlot(srt,features="Cell_Cycle_Progression",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  ggpubr::stat_compare_means(method="anova",size=8, label.y = 0.00015,label.x = 4)+xlab("")+ggtitle("L_MP1: Cell Cycle 1")+
  ylab("GS Density")
p5=VlnPlot(srt,features="CellCycle_Mitosis",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+theme(text = element_text(size = 30))+ 
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  ggpubr::stat_compare_means(method="anova",size=8, label.y = 0.00013,label.x = 4)+xlab("")+ggtitle("L_MP6: Cell Cycle 2")
p2=VlnPlot(srt,features="CellCycle_Immune_Response",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+theme(text = element_text(size = 30))+ 
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  ggpubr::stat_compare_means(method="anova",size=8, label.y = 0.00013,label.x = 4)+xlab("")+ggtitle("L_MP11: Cell Cycle Immune Response")

pdf(paste0(path_result,"validation/CellCycles_COVID19_Violin_new.pdf"),height = 6,width=30)
grid.arrange(p1,p5,p2,nrow=1,ncol=3)
dev.off()

### Cell Killing 
srt = readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD.rds"))
srtnew = CreateSeuratObject(srt@assays$RNA$counts,meta.data = srt@meta.data)
srtnew@reductions$umap = srt@reductions$umap
srt_subset = subset(srtnew,subset=cell.type.fine%in%c("CD4 T","CD4m T","CD4n T","CD8eff T","CD8m T","gd T","NK"))
Idents(srt_subset) = "cell.type.fine"
p0=DimPlot(srt_subset,group.by="cell.type.fine",label=T)+xlim(-10,2)+ylim(-6,6)
p1=FeaturePlot(srt_subset,features="IFN.gamma.induced_cytotoxicity",label=T)+xlim(-10,2)+ylim(-6,6)+labs(color = "GS Density")
p2=FeaturePlot(srt_subset,features="Cytotoxicity.Score",label=T)+xlim(-10,2)+ylim(-6,6)+labs(color = "GS Density")
pdf(paste0(path_result,"validation/CellLKilling-covid.pdf"),height = 5,width=12)
patchwork::wrap_plots(p0,p1,p2,nrow=1,ncol=3)
dev.off()

# ##################### PBMC 124 PROTEIN MARKERS #####################
# 
# # data=readRDS(paste0(path_data,"cite-seq/PBMC/PBMC_seuratobj.rds"))
# DefaultAssay(data) = "RNA"
# data <- data %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   FindNeighbors() %>%
#   RunUMAP(dims = 1:20)
# # ce <- compute.mca(object = data)
# # cells <- colnames(data)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 800)
# #
# # for (gs in 1:21){
# #   cv <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
# #   cl <- compute.cell.label(cv)
# #
# #   data@meta.data$geneset <- cv[colnames(data)]
# #   colnames(data@meta.data)[which(colnames(data@meta.data)=="geneset")] = names(gs_list_L[gs])
# #   assign(paste0("p",gs),FeaturePlot(data,
# #                                     features = names(gs_list_L[gs]),
# #                                     label=T,
# #                                     raster = T))
# # }
# # pdf(paste0(path_result,"validation/PBMC_GSDensity_L_PBMC.pdf"),height = 15, width = 35)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,
# #              nrow=3,ncol=7)
# # dev.off()
# #
# # for (gs in 1:13){
# #   cv <- run.rwr(el = el, gene_set = gs_list_M[[gs]], cells = cells)
# #   cl <- compute.cell.label(cv)
# #   data@meta.data$geneset <- cv[colnames(data)]
# #   colnames(data@meta.data)[which(colnames(data@meta.data)=="geneset")] = names(gs_list_M[gs])
# #   assign(paste0("p",gs),FeaturePlot(data,
# #                                     features = names(gs_list_M[gs]),
# #                                     label=T,
# #                                     raster = T))
# # }
# # dev.off()
# # pdf(paste0(path_result,"validation/PBMC_GSDensity_M_PBMC.pdf"),height = 20, width = 20)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,
# #              nrow=4,ncol=4)
# # dev.off()
# #
# # DefaultAssay(data) = "ADT"
# # # T cell memory
# # FeaturePlot(data,features = "CD45RO")
# # # T cell killing
# # FeaturePlot(data,features = c("CD45RA","CD8","CD3","CD30"))
# #
# # ##################### HPV PROTEIN MARKERS #####################
# #
# # # sample = c("T1","T2-CD45","T2-CD45plus","T3-CD45","T3-CD45plus")
# # # seurat = readRDS(paste0(path_data,"validation-data/HPV_HNSCC/2023/",sample[1],"/count/sample_filtered_feature_bc_matrix/seuratobj.rds"))
# # # for (s in 2:5){
# # #   seurat = merge(seurat,
# # #                  readRDS(paste0(path_data,"validation-data/HPV_HNSCC/2023/",sample[s],"/count/sample_filtered_feature_bc_matrix/seuratobj.rds")))
# # # }
# # # seurat = SCTransform(seurat,assay = "RNA",new.assay.name = "SCT_rna")
# # # seurat <- seurat %>%
# # #   FindVariableFeatures() %>%
# # #   RunPCA() %>%
# # #   FindNeighbors() %>%
# # #   RunUMAP(dims = 1:20)
# # # saveRDS(seurat,paste0(path_data,"validation-data/HPV_HNSCC/2023/citeseq-combined.rds"))
# #
# # seurat = readRDS(paste0(path_data,"validation-data/HPV_HNSCC/2023/citeseq-combined.rds"))
# #
# # DefaultAssay(seurat) = "ADT"
# #
# # pdf(paste0(path_result,"validation/HPVCiteSeq_memory.pdf"),height = 15, width = 20)
# # FeaturePlot(seurat,features=c("anti-human CD127 (IL-7R?)","anti-human CD45RO",
# #                               "anti-human CD69","anti-human CD45RA","anti-human CD27",
# #                               "anti-human CD103 (Integrin ?E)","anti-human CD49a"))
# # dev.off()
# #
# # DefaultAssay(seurat) = "RNA"
# # ce <- compute.mca(object = seurat)
# # cells <- colnames(seurat)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# #
# # for (gs in 1:21){
# #   cv <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
# #   cl <- compute.cell.label(cv)
# #
# #   seurat@meta.data$geneset <- cv[colnames(seurat)]
# #   colnames(seurat@meta.data)[which(colnames(seurat@meta.data)=="geneset")] = names(gs_list_L[gs])
# #   assign(paste0("p",gs),FeaturePlot(seurat,
# #                                     features = names(gs_list_L[gs]),
# #                                     label=T,
# #                                     raster = T))
# # }
# # pdf(paste0(path_result,"validation/HPVCiteSeq_GSDensity_L.pdf"),height = 15, width = 35)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,
# #              nrow=3,ncol=7)
# # dev.off()
# #
# # for (gs in 1:13){
# #   cv <- run.rwr(el = el, gene_set = gs_list_M[[gs]], cells = cells)
# #   cl <- compute.cell.label(cv)
# #   seurat@meta.data$geneset <- cv[colnames(seurat)]
# #   colnames(seurat@meta.data)[which(colnames(seurat@meta.data)=="geneset")] = names(gs_list_M[gs])
# #   assign(paste0("p",gs),FeaturePlot(seurat,
# #                                     features = names(gs_list_M[gs]),
# #                                     label=T,
# #                                     raster = T))
# # }
# #
# # pdf(paste0(path_result,"validation/HPVCiteSeq_GSDensity_M.pdf"),height = 20, width = 25)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,
# #              nrow=4,ncol=4)
# # dev.off()
# 
