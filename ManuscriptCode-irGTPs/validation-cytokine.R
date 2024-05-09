# library(MASS)
# library(Matrix)
# library(mclust)
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
# library(ComplexHeatmap)
# library(org.Hs.eg.db)
# library(msigdbr)
# library(fgsea)
library(forcats)
library(clusterProfiler)
library(SingleR)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

chetah.ref.singler = readRDS(paste0("/Volumes/she4/HPV_negative/data/","chetah.ref.singler.rds"))
# hpca.se = readRDS(paste0("/Volumes/she4/HPV_negative/data/","hpca.rds"))

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
gs_list = c(gs_list_L,gs_list_M)

####### Read in data

mat = readMM(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/GSE150861_matrix.mtx"))
cells = read.delim(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/GSE150861_barcodes.tsv"),header = F)
features = read.delim(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/GSE150861_features.tsv"),header = F)
mat = mat[-which(duplicated(features$V2)),]
rownames(mat) = features$V2[-which(duplicated(features$V2))]
colnames(mat) = cells$V1

metadata = read.csv(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/GSE150861_cell_annotation.csv"),row.names = 1)
srt = CreateSeuratObject(as.matrix(mat))
srt = AddMetaData(srt,metadata = metadata)

####### healthy controls
healthy = Read10X_h5(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/5k_pbmc_NGSC3_aggr_filtered_feature_bc_matrix.h5"))
healthy = CreateSeuratObject(as.matrix(healthy))
srt_healthy = merge(srt,healthy)
# srt_healthy = JoinLayers(srt_healthy)
srt_healthy$day = ifelse(srt_healthy$day%in%c("day1","day5","day7"),srt_healthy$day,"healthy")
srt_healthy <- NormalizeData(srt_healthy, verbose = F)
srt_healthy <- FindVariableFeatures(srt_healthy, verbose = F)
srt_healthy <- ScaleData(srt_healthy, verbose = F)
srt_healthy <- RunPCA(srt_healthy, features = VariableFeatures(object = srt_healthy), verbose = F)
srt_healthy <- FindNeighbors(srt_healthy, dims = 1:20, verbose = F)
srt_healthy <- FindClusters(srt_healthy, resolution = 1., verbose = F)
srt_healthy <- RunUMAP(srt_healthy, dims = 1:20,seed.use = 42, verbose = F)
saveRDS(srt_healthy,paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/srt_healthy.rds"))

srt_healthy = readRDS(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/srt_healthy.rds"))

# FeaturePlot(srt_healthy,features=c("CD14","CD163","LYZ","CD68"),raster = F,label=T)
FeaturePlot(srt_healthy,features=c("CD14","LYZ","CCL4","IL1B",
                                   "TNF","CCL3","CXCL2","IL10"),raster = F,label=T)
type = SingleR(test = as.SingleCellExperiment(srt_healthy), ref = chetah.ref.singler[["data"]],
               assay.type.test=1,labels = chetah.ref.singler$types)
srt_healthy = AddMetaData(srt_healthy,metadata=type$labels,col.name='celltype_snglr_chetah')
DimPlot(srt_healthy,group.by = c("celltype_snglr_chetah","day"),label = T)

srt_subset = subset(srt_healthy,subset=seurat_clusters%in%c("3","8","10","18","23"))
DimPlot(srt_subset,group.by = c("celltype_snglr_chetah","day","seurat_clusters"),
        label = T)

srt_subset$compare = ifelse(srt_subset$seurat_clusters%in%c("23","18"),
                            "ICS","Other_Macrophage")
Idents(srt_subset) = "compare"
# saveRDS(srt_subset,paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/srt_healthy_subset.rds"))

# srt_subset = subset(srt,cells=colnames(srt)[which(srt$celltype_snglr_chetah=="Macrophage")])
# srt_subset <- RunPCA(srt_subset, features = VariableFeatures(object = srt_subset), verbose = F)
# srt_subset <- FindNeighbors(srt_subset, dims = 1:20, verbose = F)
# srt_subset <- FindClusters(srt_subset, resolution = 1., verbose = F)
# srt_subset <- RunUMAP(srt_subset, dims = 1:20,seed.use = 42, verbose = F)

# srt_subset = readRDS(paste0(path_data,"validation-data/cytokine/COVID19-cytokinestorm/srt_healthy_subset.rds"))
# Idents(srt_subset) = "day"
# srt_subset1 = subset(srt_subset,downsample=900)
# DimPlot(srt_healthy,group.by = "seurat_clusters",label = T)
# cells_AUC <- AUCell_run(srt_healthy@assays$RNA$counts, geneSets=list(geneset = gs_list_L[[18]]))
# srt_healthy$score = cells_AUC@assays@data@listData[["AUC"]]
# FeaturePlot(srt_healthy,features="score")
# VlnPlot(subset(srt_healthy,subset=celltype_snglr_chetah=="Macrophage"),
#         features="score",group.by = "seurat_clusters",pt.size=0)

srt_subset = subset(srt_healthy,subset=celltype_snglr_chetah%in%c("Macrophage"))
srt_subset = CreateSeuratObject(srt_subset@assays$RNA$counts,meta.data = srt_subset@meta.data)
Idents(srt_subset) = "day"
# srt_subset1 = subset(srt_subset,downsample=2000)
markers = FindMarkers(srt_subset,ident.1="day1",ident.2="day5",logfc.threshold = 0)
gene_marker_subset = markers
ranks <- gene_marker_subset$avg_log2FC
names(ranks) <- rownames(gene_marker_subset)
pvalue = fgsea(pathways = gs_list, stats = ranks)

p0=plotEnrichment(gs_list[[6]],ranks)+ labs(title="Chemokine_Mediated_T_activation")
p0$layers[[5]]$aes_params$colour <- 'purple'
p0 = p0+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Chemokine_Mediated_T_activation"),6],3)," Adjusted p-value: ",
                                 round(pvalue[which(pvalue$pathway == "Chemokine_Mediated_T_activation"),3],3)))+
    theme(text = element_text(size = 8,face="bold"))
  
p1=plotEnrichment(gs_list[[3]],ranks)+ labs(title="Cytokine_Production")
p1$layers[[5]]$aes_params$colour <- 'purple'
p1 = p1+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Cytokine_Production"),6],3)," Adjusted p-value: ",
                                 round(pvalue[which(pvalue$pathway == "Cytokine_Production"),3],3)))+
    theme(text = element_text(size = 8,face="bold"))
  
p2=plotEnrichment(gs_list[[9]],ranks)+ labs(title="Chemokine_activity")
p2$layers[[5]]$aes_params$colour <- 'purple'
p2 = p2+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Chemokine_activity"),6],3)," Adjusted p-value: ",
                                   round(pvalue[which(pvalue$pathway == "Chemokine_activity"),3],3)))+
      theme(text = element_text(size = 8,face="bold"))
    
p3=plotEnrichment(gs_list[[27]],ranks)+ labs(title="Cytokine_Receptor_Signaling")
p3$layers[[5]]$aes_params$colour <- 'purple'
p3 = p3+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Cytokine_Receptor_Signaling"),6],3)," Adjusted p-value: ",
                                     round(pvalue[which(pvalue$pathway == "Cytokine_Receptor_Signaling"),3],3)))+
        theme(text = element_text(size = 8,face="bold"))
pdf(paste0(path_result,"validation/cytokinestorm_barcode.pdf"),height=8,width = 5)
patchwork::wrap_plots(p1,p2,p3,nrow=3)
dev.off()
pdf(paste0(path_result,"validation/cytokinestorm_barcode_horizontal.pdf"),height=3,width = 12)
patchwork::wrap_plots(p0,p1,p2,p3,nrow=1)
dev.off()
      
# Cytokine
# gene_marker_subset = markers[order(markers$avg_log2FC,decreasing = T),]
# gene_list = gene_marker_subset$avg_log2FC
# names(gene_list) = rownames(gene_marker_subset)
# gsea = GSEA(gene_list,TERM2GENE = TERM2GENE,pvalueCutoff=1,minGSSize=0)
# gsea = gsea@result
# gsea = gsea[order(gsea$NES,decreasing = T),]
# gsea$Description2 <- factor(gsea$Description, levels = rev(rownames(gsea)))
# gsea$psig = ifelse(gsea$p.adjust<=0.05,"Sig","Non-sig")
# pdf("/rsrch4/home/bcb/she4/hallmark/results/validation/cytokine_covidstorm_alldotplot.pdf",
#     height=6,width = 9)
# ggplot(gsea, aes(x = NES, y = Description2)) +
#   geom_point(aes(size = NES, color = psig)) +
#   theme_bw(base_size = 14) +
#   scale_colour_manual(values=c("grey","red")) +
#   ylab(NULL) +xlab("Normalized Enrichment Score")+
#   ggtitle("Severe COVID-19 vs. Healthy")+theme(text=element_text(size=16,face = "bold"))
# dev.off()

# annotation_M = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))
# TERM2GENE_M = NULL
# file_names_MP <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
# for (i in 1:length(file_names_MP)){
#   TERM2GENE_M = rbind(TERM2GENE_M,cbind(rep(names(annotation_M)[i],50),read.csv(file_names_MP[i])[,2]))
# }
# TERM2GENE_M = as.data.frame(TERM2GENE_M)
# names(TERM2GENE_M) = c("term",'gene')
# annotation_L = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
# names(annotation_L)[15] = "Cell_Adhesion_1_L"
# TERM2GENE_L = NULL
# file_names_MP <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
# for (i in 1:length(file_names_MP)){
#   TERM2GENE_L = rbind(TERM2GENE_L,cbind(rep(names(annotation_L)[i],50),read.csv(file_names_MP[i])[,2]))
# }
# TERM2GENE_L = as.data.frame(TERM2GENE_L)
# names(TERM2GENE_L) = c("term",'gene')
# TERM2GENE = as.data.frame(rbind(TERM2GENE_M,TERM2GENE_L))
# 
# for (c in 1:length(unique(markers$cluster))){
#   gene_marker_subset = markers[which(markers$cluster==unique(markers$cluster)[c]),]
#   if(nrow(gene_marker_subset)==0){next}
#   gene_marker_subset = gene_marker_subset[order(gene_marker_subset$avg_log2FC,decreasing = T),]
#   gene_list = gene_marker_subset$avg_log2FC
#   names(gene_list) = gene_marker_subset$gene
#   gsea = GSEA(gene_list,TERM2GENE = TERM2GENE,pvalueCutoff=1,minGSSize=0)
#   write.csv(as.data.frame(gsea@result),paste0(path_result,"validation/cytokine_gsea_COVIDstorm_",unique(markers$cluster)[c],".csv"))
# }
# 
# # Cytokine
# gsea1 = read.csv("/Volumes/she4-1/hallmark/results/validation/cytokine_gsea_COVIDstorm_day1.csv",row.names = 1)
# gsea1 = gsea1[order(gsea1$NES,decreasing = T),]
# gsea2 = read.csv("/Volumes/she4-1/hallmark/results/validation/cytokine_gsea_COVIDstorm_day7.csv",row.names = 1)
# gsea2 = gsea2[rownames(gsea1),]
# gsea = as.data.frame(rbind(gsea1,gsea2))
# gsea$COVID_severity = as.factor(c(rep("Severe COVID-19 patients",28),rep("Tocilizumab-treated",28)))
# gsea$Description2 <- factor(gsea$Description, levels = rev(rownames(gsea1)))
# 
# pdf("/rsrch4/home/bcb/she4/hallmark/results/validation/cytokine_covidstorm.pdf",
#     height=6,width = 9)
# ggplot(gsea[which(gsea$COVID_severity=="Severe COVID-19 patients"),], aes(x = NES, y = Description2)) +
#   geom_point(aes(size = NES, color = p.adjust)) +
#   theme_bw(base_size = 14) +
#   scale_colour_gradient(limits=c(0, 0.10), low="red") +
#   ylab(NULL) +xlab("Normalized Enrichment Score")+facet_grid(~COVID_severity)+
#   ggtitle("Severe COVID-19 vs. Tocilizumab-treated")+theme(text=element_text(size=16,face = "bold"))
# dev.off()
# 
# markers = read.csv(paste0(path_result,"cytokinestorm_gsea_COVIDstorm_healthy.csv"),row.names = 1)
# gene_marker_subset = markers
# ranks <- gene_marker_subset$avg_log2FC
# names(ranks) <- rownames(gene_marker_subset)
# gs_list = c(gs_list_M,gs_list_L)
# pvalue = fgsea(pathways = gs_list, stats = ranks)
# 
# p1=plotEnrichment(gs_list[[3]],ranks)+ labs(title="Cytokine_Production")
# p1$layers[[5]]$aes_params$colour <- 'purple'
# p1 = p1+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Cytokine_Production"),6],3)," Adjusted p-value: ",
#                                  round(pvalue[which(pvalue$pathway == "Cytokine_Production"),3],3)))+
#     theme(text = element_text(size = 10,face="bold"))
#   
# p2=plotEnrichment(gs_list[[9]],ranks)+ labs(title="Chemokine_activity")
# p2$layers[[5]]$aes_params$colour <- 'purple'
# p2 = p2+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Chemokine_activity"),6],3)," Adjusted p-value: ",
#                                    round(pvalue[which(pvalue$pathway == "Chemokine_activity"),3],3)))+
#       theme(text = element_text(size = 10,face="bold"))
#     
# p3=plotEnrichment(gs_list[[27]],ranks)+ labs(title="Cytokine_Receptor_Signaling")
# p3$layers[[5]]$aes_params$colour <- 'purple'
# p3 = p3+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Cytokine_Receptor_Signaling"),6],3)," Adjusted p-value: ",
#                                      round(pvalue[which(pvalue$pathway == "Cytokine_Receptor_Signaling"),3],3)))+
#         theme(text = element_text(size = 10,face="bold"))
# pdf(paste0(path_result,"validation/cytokinestorm_barcode.pdf"),height=8,width = 5)
# patchwork::wrap_plots(p1,p2,p3,nrow=3)
# dev.off()
# 
