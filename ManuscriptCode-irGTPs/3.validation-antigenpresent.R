library(MASS)
library(mclust)
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
library(ggplot2)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"
  
######### Calculate pTMB : number of single copy mutation

# copy = read.delim(paste0(path_data,"validation-data/antigenpresent/TCGA_mastercalls.abs_segtabs.fixed.txt"))
# mutation = read.delim(paste0(path_data,"validation-data/antigenpresent/mc3.v0.2.8.PUBLIC.maf"))
# 
# pTMB = NULL
# copy_subset = copy[which(copy$Modal_Total_CN!="0"),]
# for (i in 1:nrow(copy_subset)){
#   bool = mutation[,"Start_Position"]>=copy_subset[i,"Start"]&mutation[,"End_Position"]<=copy_subset[i,"End"]&copy_subset$Chromosome[i]==mutation$Chromosome
#   if(sum(bool)==0){
#     pTMB = c(pTMB,0)
#   }else{
#     pTMB = c(pTMB,length(which(bool=="TRUE")))
#   }
# }
# 
# dat = as.data.frame(cbind(copy_subset[,"Sample"],pTMB))
# dat_bulk = aggregate(x=as.numeric(dat$pTMB),by=list(dat$V1),FUN="sum")
# dat_bulk$pTMB_group = ifelse(dat_bulk$x<mean(dat_bulk$x),1,2)
# 
# # gm = Mclust(dat_bulk$x, G=2, model="V")
# # table(gm$classification)
# # dat_bulk$pTMB_group = gm$classification
# 
# ######### Create Seurat Object 
# 
# count = read.delim(paste0(path_data,"validation-data/antigenpresent/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"))
# rownames(count) = count[,1]
# count = count[,-1]
# colnames(count) = substr(colnames(count),1,12)
# colnames(count) = gsub("\\.", "-", colnames(count))
# dat_bulk$Group.1 = substr(dat_bulk$Group.1,1,12)
# 
# # Metadata
# names(dat_bulk) = c("Sample","pTMB_count","pTMB_group")
# dat_bulk = dat_bulk[!duplicated(dat_bulk$Sample),]
# rownames(dat_bulk) = dat_bulk$Sample
# 
# common_sample = intersect(dat_bulk$Sample,colnames(count))
# 
# # Expression Matrix
# mat = count[,common_sample]
# row_names = sub("\\|.*", "", rownames(mat))
# mat = mat[-which(row_names=="?"),]
# mat = mat[-which(duplicated(row_names[-which(row_names=="?")])=="TRUE"),]
# rownames(mat) = row_names[-which(row_names=="?")][-which(duplicated(row_names[-which(row_names=="?")])=="TRUE")]
# mat[is.na(mat)] = 0
# 
# # Create Seurat Object
# srt = CreateSeuratObject(counts = mat, project = "pTMB")
# srt <- AddMetaData(object = srt,metadata = dat_bulk)
# Idents(srt) = "pTMB_group"
# saveRDS(srt, paste0(path_data,"validation-data/TCGA/antigen_pTMB.rds"))

######### DEGs and GSEA

srt = readRDS(paste0(path_data,"/validation-data/TCGA/antigen_pTMB.rds"))
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
srt <- ScaleData(srt, features = rownames(srt))
markers = FindAllMarkers(srt,logfc.threshold=0,return.thresh = 1,min.pct = 0)
write.csv(markers,paste0(path_result,"antigenpresent_pTMB_markers.csv"))

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

for (c in 1:length(unique(markers$cluster))){
  gene_marker_subset = markers[which(markers$cluster==unique(markers$cluster)[c]),]
  if(nrow(gene_marker_subset)==0){next}
  gene_marker_subset = gene_marker_subset[order(gene_marker_subset$avg_log2FC,decreasing = T),]
  gene_list = gene_marker_subset$avg_log2FC
  names(gene_list) = gene_marker_subset$gene
  gsea = GSEA(gene_list,TERM2GENE = TERM2GENE,pvalueCutoff=1,minGSSize=0)
  write.csv(as.data.frame(gsea@result),paste0(path_result,"validation/antigenpresent_gsea_M_pTMB_",c,"scale.csv"))
}
for (c in 1:length(unique(markers$cluster))){
  gene_marker_subset = markers[which(markers$cluster==unique(markers$cluster)[c]),]
  if(nrow(gene_marker_subset)==0){next}
  gene_marker_subset = gene_marker_subset[order(gene_marker_subset$avg_log2FC,decreasing = T),]
  gene_list = gene_marker_subset$avg_log2FC
  names(gene_list) = gene_marker_subset$gene
  gsea = GSEA(gene_list,TERM2GENE = TERM2GENE_L,pvalueCutoff=1,minGSSize=0)
  write.csv(as.data.frame(gsea@result),paste0(path_result,"validation/antigenpresent_gsea_L_pTMB_",c,"scale.csv"))
}

markers = read.csv(paste0(path_result,"antigenpresent_pTMB_markers.csv"),row.names = 1)
gene_marker_subset = markers[which(markers$cluster==1),]
ranks <- gene_marker_subset$avg_log2FC
names(ranks) <- gene_marker_subset$gene

pvalue = fgsea(pathways = gs_list_L, stats = ranks)

p1=plotEnrichment(gs_list_L[[2]],ranks)+ labs(title="TCR Anchoring")
p1$layers[[5]]$aes_params$colour <- 'purple'
p1 = p1+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "TCR_Anchoring"),6],3)," Adjusted p-value: ",
                          round(pvalue[which(pvalue$pathway == "TCR_Anchoring"),3],3)))+
  theme(text = element_text(size = 10))
  
p2=plotEnrichment(gs_list_L[[3]],ranks)+ labs(title="L_MP3: Histone_associated_lipid_antigen_presentation")
p2$layers[[5]]$aes_params$colour <- 'purple'
p2 = p2+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Histone_associated_lipid_antigen_presentation"),6],3)," Adjusted p-value: ",
                           round(pvalue[which(pvalue$pathway == "Histone_associated_lipid_antigen_presentation"),3],3)))+
  theme(text = element_text(size = 10,face="bold"))

p3=plotEnrichment(gs_list_L[[5]],ranks)+ labs(title="L_MP5: Antigen Presentation")
p3$layers[[5]]$aes_params$colour <- 'purple'
p3 = p3+labs(subtitle = paste0("NES: ",round(pvalue[which(pvalue$pathway == "Antigen_Presentation"),6],3)," Adjusted p-value: ",
                            round(pvalue[which(pvalue$pathway == "Antigen_Presentation"),3],3)))+
  theme(text = element_text(size = 10,face="bold"))

pdf(paste0(path_result,"validation/antigen_pTMB_barcode.pdf"),height=7,width = 6)
grid.arrange(p2,p3,nrow=2)
dev.off()
pdf(paste0(path_result,"validation/antigen_pTMB_barcode_LMP3.pdf"),height=3,width = 4)
p2
dev.off()
pdf(paste0(path_result,"validation/antigen_pTMB_barcode_LMP5.pdf"),height=3,width = 4)
p3
dev.off()


markers = read.csv(paste0(path_result,"antigenpresent_pTMB_markers.csv"),row.names = 1)
for (c in 1:length(unique(markers$cluster))){
  gene_marker_subset = markers[which(markers$cluster==unique(markers$cluster)[c]),]
  if(nrow(gene_marker_subset)==0){next}
  gene_marker_subset = gene_marker_subset[order(gene_marker_subset$avg_log2FC,decreasing = T),]
  gene_list = gene_marker_subset$avg_log2FC
  names(gene_list) = gene_marker_subset$gene
  gsea = GSEA(gene_list,TERM2GENE = TERM2GENE_L,pvalueCutoff=1,minGSSize=0)
  # write.csv(as.data.frame(gsea@result),paste0(path_result,"validation/celladhesion_gsea_M_TCGA_",c,"_scaled.csv"))
}
gsea = gsea[order(gsea$NES,decreasing = T),]
gsea$Description <- factor(gsea$Description, levels = rev(rownames(gsea)))
pdf(paste0(path_result,"validation/antigen_pTMB_all.pdf"),height=6,width = 9)
ggplot(gsea, aes(x = NES, y = Description)) +
  geom_point(aes(size = NES, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 1), low="brown2",high="cornflowerblue") +
  ylab(NULL) +xlab("Normalized Enrichment Score")+
  ggtitle("Enrichment for high pTMB samples")+
  theme(text=element_text(size=13,face = "bold"))
dev.off()
