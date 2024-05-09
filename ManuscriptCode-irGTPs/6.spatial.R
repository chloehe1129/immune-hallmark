library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(rhdf5)
library(CellTrek)
library(semla)

options(stringsAsFactors = F)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_data = "/rsrch4/home/bcb/she4/hallmark/results/"

gs_list_L = list()
file_names_L <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names_L)){
  gs_list_L[[i]] = read.csv(file_names_L[i],row.names = 1)[,1]
}
annotation_L = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
names(gs_list_L) = paste0("L_",names(annotation_L))

gs_list_M = list()
file_names_M <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names_M)){
  gs_list_M[[i]] = read.csv(file_names_M[i],row.names = 1)[,1]
}
annotation_M = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))
names(gs_list_M) = paste0("M_",names(annotation_M))

exhaustion=c('PDCD1','TIGIT','HAVCR2','CTLA4','LAG3') #from zeming zhang
effector=c('NKG7','PRF1','GZMA','GZMB','GZMH','GNLY',"IFNG")
exh = read.delim(paste0(path_data,"exh.txt"),header = F)
exh = exh$V1
naive=c('CCR7','LEF1','SELL','TCF7') #from CXCR6/CXCL16 paper
tissue_memory.full=c('ITGAE','ITGA1','ZNF683','CD69','IL2','IL10',"CXCL6","CXCL13",
                     'CA10','CXCR6','CXCL13','DUSP6','KCNK5','RGS1','CRTAM','PDCD1')
TEM = c("CX3CR1","GZMB","GNLY","GZMK","GZMM","CD44")
Tprolif = c("MKI67","TUBB","STMN1")
CD8 = c("CD8A","CD69")
Treg = c("FOXP3","IL2RA")
M1 = c("CD80","CD86","CD68","IL1R1","TLR2","TLR4","SOCS3","TNF")
M2 = c("CD163","MRC1","CD200R1","TGM2","IL1R2","TGFB1","IL10")
genesetlist=list(exh,exhaustion,effector,naive,tissue_memory.full,TEM,Tprolif,CD8,Treg,M1,M2)
genesetlist_names=c("Exhaustion",'exhaustion_short','Cytotoxicity',
                    'naive','Tissue.Resident.Memory',"Tissue.Effector.Memory","Tprolif","CD8","Treg","M1",
                    "M2")
names(genesetlist) = genesetlist_names

# Hallmark
gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
H_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(H_list) = unique(gene_sets$gs_name)

# KEGG
kegg_sets = read.csv(paste0("/Volumes/she4/hallmark/data/","kegg_immune_sets.csv"))[,-1]
KEGG_list = list()
for (i in 1:length(unique(kegg_sets$name))){
  KEGG_list[[i]] = unique(kegg_sets[which(kegg_sets$name==unique(kegg_sets$name)[i]),"symbol"])
}
names(KEGG_list) = unique(kegg_sets$name)

gs_list = c(gs_list_L,gs_list_M,genesetlist,H_list,KEGG_list)

TLS = c("FDCSP","CR2","CXCL13","LTF","CD52","MS4A1","CCL19","
          LINC00926","LTB","CORO1A","CD79B","TXNIP","CD19","LIMD2","CD37","ARHGAP45","BLK","TMC8","
          CCL21","PTPN6","ATP2A3","IGHM","SPIB","TMSB4X","CXCR4","NCF1","CD79A","ARHGAP9","DEF6","EVL","TBC1D10C","RASAL3","
          INPP5D","RNASET2","RASGRP2","TNFRSF13C","RAC2","CD22","ARHGEF1","AC103591.3","TRAF3IP3","HLA-DQB1","CD53","ARHGAP4","TRBC2","POU2AF1","TRAF5","OGA","FCRL3","HLA-DQA1")
####################### Supplementary to Figure 6 #######################
brast_ST = Load10X_Spatial(paste0(path_data,"spatial/breast/"))
# domain = read.csv(paste0(path_data,"spatial/breast/final.domain.csv"),row.names = 1)
# brast_ST = AddMetaData(brast_ST,metadata=domain)
for (j in 1:length(gs_list)){
  marker_list <- gs_list[[j]][gs_list[[j]] %in% rownames(brast_ST@assays$Spatial$counts)]
  print (j); length(marker_list)
  name_new = names(gs_list)[j]
  brast_ST <- AddModuleScore(brast_ST, features = list(marker_list), nbin = 12, name = name_new)
}
TLS = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/10X_pathologist/breast/TLS.csv",
               row.names = 1)
TLS$TLS = ifelse(TLS$TLS=="TLS","TLS","Not TLS")
brast_ST = AddMetaData(brast_ST,metadata = TLS)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_TLS_breast.pdf",
    height=4,width=4)
VlnPlot(brast_ST,features = "L_TCR_Anchoring1",group.by = "TLS",pt.size=0,
        cols = c("hotpink","cornflowerblue"))+
  ggpubr::stat_compare_means(label="p.signif")+ylim(-2.5,5)
dev.off()

malignancy = read.delim(paste0(path_data,"spatial/breast/BC_copykat_prediction.txt"))
rownames(malignancy) = malignancy$cell.names
brast_ST = AddMetaData(brast_ST,metadata = malignancy)

pseudotime = read.csv(paste0(path_data,"spatial/breast/pseudotime.csv"),row.names = 1)
pseudotime[is.na(pseudotime)] = 0
brast_ST = AddMetaData(brast_ST,metadata = pseudotime)

####################### Correlation between pseudotime and MP activities #######################
# cluster 6 - myeloid
# cluster 7 - lymphoid
cor = cor(brast_ST@meta.data[,c(4:22,125)])
p = ggcorrplot::cor_pmat(brast_ST@meta.data[,c(4:22,125)])
cor = cor[,20]
p = p[,20]
p = p.adjust(p,method="fdr")
p1=Heatmap(as.matrix(cor),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(paste0(as.matrix(round(cor,2))[i, j],ifelse(as.matrix(p)[i, j]<0.05,"*","")), x, y, gp = gpar(fontsize = 10))
           },show_row_dend = F,name="Correlation",row_names_max_width = unit(8,"cm"))

cor = cor(brast_ST@meta.data[,c(23:31,127)])
p = ggcorrplot::cor_pmat(brast_ST@meta.data[,c(23:31,127)])
cor = cor[,10]
p = p[,10]
p = p.adjust(p,method="fdr")
p2=Heatmap(as.matrix(cor),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(paste0(as.matrix(round(cor,2))[i, j],ifelse(as.matrix(p)[i, j]<0.05,"*","")), x, y, gp = gpar(fontsize = 10))
           },show_row_dend = F,name="Correlation",row_names_max_width = unit(8,"cm"))

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_L_MPvsPseudotime_breast.pdf",
    height=8,width=6)
p1
dev.off()
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_M_MPvsPseudotime_breast.pdf",
    height=4,width=5)
p2
dev.off()

######## Incorporate CytoSpace celltype deconvolution ########
cell_type_assignments_by_spot = read.csv("/Volumes/she4/hallmark/results/spatial/cytospace_results/cell_type_assignments_by_spot.csv",row.names = 1)
cell_type_assignments_by_spot$T_percent = (cell_type_assignments_by_spot$T.cells.CD8+cell_type_assignments_by_spot$T.cells.CD4)/cell_type_assignments_by_spot$Total.cells
cell_type_assignments_by_spot$B_percent = (cell_type_assignments_by_spot$B.cells)/cell_type_assignments_by_spot$Total.cells
cell_type_assignments_by_spot$Myeloid_percent = (cell_type_assignments_by_spot$Monocytes.and.Macrophages)/cell_type_assignments_by_spot$Total.cells
cell_type_assignments_by_spot$NK_percent = (cell_type_assignments_by_spot$NK.cells)/cell_type_assignments_by_spot$Total.cells
cell_type_assignments_by_spot$DC_percent = (cell_type_assignments_by_spot$Dendritic.cells)/cell_type_assignments_by_spot$Total.cells
cell_type_assignments_by_spot$barcode = unlist(strsplit(rownames(cell_type_assignments_by_spot),"__"))[seq(1,nrow(cell_type_assignments_by_spot)*2,2)]
cell_type_assignments_by_spot$barcode = gsub("\\.","-",cell_type_assignments_by_spot$barcode)
rownames(cell_type_assignments_by_spot) = cell_type_assignments_by_spot$barcode
brast_ST = AddMetaData(brast_ST,metadata = cell_type_assignments_by_spot)

########## Tumor Boundaries and gene set relationship ##########
brast_ST = UpdateSeuratForSemla(object=brast_ST,verbose = TRUE)
brast_ST <- SetIdent(brast_ST, value = "cell_names")
for(s in 1:length(unique(brast_ST$cell.names))){ # zero neighbors: 1715
  if(s == 1715)next
  brast_ST <- RegionNeighbors(brast_ST, column_name = "cell.names",
                              column_labels = unname(unique(brast_ST$cell.names)[s]), 
                              verbose = TRUE)
  print(s)
}
spot = names(brast_ST@meta.data)[142:ncol(brast_ST@meta.data)]
neighbor_res = matrix(0,length(spot),1)
rownames(neighbor_res) = spot
for(s in spot){
  neighbors = rownames(brast_ST@meta.data)[which(!is.na(brast_ST@meta.data[,s]))]
  malignant_percent = sum(malignancy[neighbors,2]=="aneuploid")/length(neighbors)
  neighbor_res[s,] = malignant_percent
}
rownames(neighbor_res) = substr(rownames(neighbor_res),7,nchar(rownames(neighbor_res)))
names(neighbor_res) = "Malignant_percent"
brast_ST = AddMetaData(brast_ST,metadata = as.data.frame(neighbor_res))

# brast_ST_T = subset(brast_ST,subset=T_percent>0.8&copykat.pred=="diploid")
brast_ST_T = subset(brast_ST,subset=Myeloid_percent>0.6&copykat.pred=="diploid")
tmp = brast_ST_T@meta.data[,c(4:92,which(names(brast_ST_T@meta.data)=="V1"))]
# tmp$M1_M2_ratio = tmp$M11/tmp$M21
tmp = tmp[,c(5,14,15,29,31,32,33,34,90)] # For T spots 5,14,15
# tmp = tmp[,c(20:28,83:84,91,90)] # For Myeloid spots
mat = NULL
for(i in 1:8){
  mat = rbind(mat,aggregate(tmp[,i],by=list(tmp$V1),FUN=mean)[,2])
}
mat1 = apply(mat,1,scale)
colnames(mat1) = names(tmp)[1:8]
rownames(mat1) = as.character(round(sort(unique(tmp$V1)),2))
anno = matrix(0,8,2)
options(scipen = 0)
for(i in 1:8){
  anno[i,1]=round(cor.test(tmp$V1,tmp[,1:8][,i])$estimate,2)
  anno[i,2]=format.pval(cor.test(tmp$V1,tmp[,1:8][,i])$p.value)
}
rownames(anno) = names(tmp)[1:8]
anno[,1] = as.numeric(anno[,1])
anno[,2] = as.numeric(ifelse(anno[,2]=="<2e-16","0",anno[,2]))
anno = anno[order(anno[,1],decreasing = T),]
row_ha = circlize::colorRamp2(seq(min(as.numeric(anno[,1])),max(as.numeric(anno[,1])),0.1),hcl_palette = "RdBu",reverse = T)
pval = symnum(as.numeric(anno[,2]), corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "..", "."))
row_ha = rowAnnotation(Correlation = as.numeric(anno[,1]),col=list(Correlation=row_ha),
                       Significance = anno_text(paste0(unname(anno[,1]),pval)),show_annotation_name=F)
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "violetred3"))

col_ha = HeatmapAnnotation(Tumor.Proximity = anno_lines(as.numeric(rownames(mat1))))

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_boundary_vs_GS_heatmap_Myeloid.pdf",
    height=4,width=9)
Heatmap(t(mat1)[rownames(anno),],name = "Z-norm gene set score",col=col_fun,
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = col_ha,
        show_column_dend = F,show_row_dend = F,
        row_names_max_width = unit(10, "cm"),
        left_annotation = row_ha,
        # km=2,
        border = T,
        row_split = ifelse(as.numeric(anno[,1])>=0,"positive","negative"),
        column_title = "Neighboring Malignant %",column_title_side = "bottom")
dev.off()

# tmp = brast_ST_T@meta.data[-which(brast_ST_T@meta.data$V1==1|brast_ST_T@meta.data$copykat.pred=="aneuploid"),c(4:90,which(names(brast_ST_T@meta.data)=="V1"))]

# signatures = c("HALLMARK_ALLOGRAFT_REJECTION1","HALLMARK_INFLAMMATORY_RESPONSE1",
#                "HALLMARK_INTERFERON_GAMMA_RESPONSE1",
#                "naive1","Treg1","Tissue.Resident.Memory1",
#                "Cytotoxicity1","Tissue.Effector.Memory1","Exhaustion1")
# signatures = c("HALLMARK_ALLOGRAFT_REJECTION1","HALLMARK_INFLAMMATORY_RESPONSE1",
#                "HALLMARK_INTERFERON_GAMMA_RESPONSE1","HALLMARK_INTERFERON_ALPHA_RESPONSE1",
#                "M11","M21")
# temp = tmp[,c("V1",signatures)]
# temp$V1 = round(temp$V1,2)
# temp$V1 = factor(temp$V1,sort(unique(temp$V1)))
# temp = aggregate(temp[,signatures],by=list(temp$V1),FUN=mean)
# temp_long = melt(setDT(temp),id.var="Group.1")
# temp_long$value = round(temp_long$value,2)
# temp_long$variable = factor(temp_long$variable,signatures)
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_Infiltration_vs_Distance_Myeloid.pdf",
#     height=4,width=11)
# ggplot(temp_long,aes(x=Group.1,y=value))+
#   geom_histogram(stat = "identity",position = position_dodge())+
#   facet_wrap(~variable,scale="free_y",ncol=3)+
#   xlab("Neighboring Malignant %")+ylab("Z-norm Cell Type Score")+
#   geom_hline(yintercept=0, linetype="dashed", color = "brown3", size=0.6)+theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# dev.off()

# mat = t(as.matrix(aggregate(tmp[,c("V1","CD81","effector1")],by=list(tmp$V1),FUN=mean)[,c(3,4)]))
# mat = t(apply(mat,1,scale))
# colnames(mat) = as.factor(round(sort(unique(tmp$V1)),2))
# pdf("~/desktop/infiltration_down_heatmap.pdf",height=1,width=10)
# Heatmap(mat,name = "Z-norm score",col=col_fun,show_column_dend = F,show_row_dend = F)
# dev.off()
# mat = t(as.matrix(aggregate(tmp[,c("V1","HALLMARK_G2M_CHECKPOINT1")],by=list(tmp$V1),FUN=mean)[,c(3)]))
# mat = t(apply(mat,1,scale))
# colnames(mat) = as.factor(round(sort(unique(tmp$V1)),2))
# rownames(mat) = "G2M"
# pdf("~/desktop/infiltration_up_G2M.pdf",height=1,width=10)
# Heatmap(mat,name = "Z-norm score",col=col_fun,show_column_dend = F,show_row_dend = F,
#         row_names_max_width = unit(10, "cm"))
# dev.off()
mat = NULL
for(i in 1:19){
  mat = rbind(mat,aggregate(tmp[,i],by=list(tmp$V1),FUN=mean)[,2])
}
mat1 = apply(mat,1,scale)
colnames(mat1) = names(tmp)[1:19]
rownames(mat1) = as.character(round(sort(unique(tmp$V1)),2))
anno = matrix(0,19,2)
for(i in 1:19){
  anno[i,1]=round(cor.test(tmp$V1,tmp[,1:19][,i])$estimate,2)
  anno[i,2]=format.pval(cor.test(tmp$V1,tmp[,1:19][,i])$p.value)
}
rownames(anno) = names(tmp)[1:19]
anno[,1] = as.numeric(anno[,1])
anno[,2] = as.numeric(anno[,2])
anno = anno[order(anno[,1],decreasing = T),]
row_ha = circlize::colorRamp2(seq(min(as.numeric(anno[,1])),max(as.numeric(anno[,1])),0.1),hcl_palette = "RdBu",reverse = T)
pval = symnum(as.numeric(anno[,2]), corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
row_ha = rowAnnotation(Correlation = as.numeric(anno[,1]),col=list(Correlation=row_ha),Significance = anno_text(pval))
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "violetred3"))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_boundary_vs_GS_heatmap_Myeloid.pdf",
    height=4,width=12)
Heatmap(t(mat1)[rownames(anno),],name = "Z-norm MP score",col=col_fun,
        cluster_rows = F,
        cluster_columns = F,
        show_column_dend = F,show_row_dend = F,
        row_names_max_width = unit(10, "cm"),
        left_annotation = row_ha,
        # km=2,
        border = T,
        row_split = ifelse(as.numeric(anno[,1])>0,"positive","negative"),
        column_title = "Neighboring Malignant %",column_title_side = "bottom")
dev.off()
library(GGally)
lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping)+ 
    geom_point(alpha = .25) + 
    geom_smooth(method = "lm", formula = y ~ x, 
                fill = "blue", color = "red", size = 0.5)
}
# Myeloid
ggpairs(tmp, columns = c(20:22,24:26,28,38,64:66),lower = list(continuous = wrap(lowerfun)))+
  theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))
ggpairs(tmp, columns = c(23,27,38,39,64:66),lower = list(continuous = wrap(lowerfun)))+
  theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))

corr = cor(tmp[,c(20:22,24:26,28,38,64:66)])
p.mat = cor_pmat(tmp[,c(20:22,24:26,28,38,64:66)])
p1=ggcorrplot(corr, lab = TRUE,colors = c("cornflowerblue", "white", "salmon"))+ggtitle("Negative correlation with tumor proximity")

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_correlation_Myeloid_breast.pdf",
    height=9,width=10)
p1
dev.off()


# T
corr = cor(tmp[,c(5,14,15,29,31)])
p.mat = cor_pmat(tmp[,c(5,14,15,29,31)])
# ggcorrplot(corr, p.mat = p.mat,insig="blank",lab = TRUE)
p1=ggcorrplot(corr, lab = TRUE,colors = c("cornflowerblue", "white", "salmon"))+ggtitle("Negative correlation with tumor proximity")

corr = cor(tmp[,c(4,17,2,32,33,37)])
p.mat = cor_pmat(tmp[,c(4,17,2,32,33,37)])
# ggcorrplot(corr, p.mat = p.mat,insig="blank",lab = TRUE)
p2=ggcorrplot(corr, lab = TRUE,colors = c("cornflowerblue", "white", "salmon"))+ggtitle("Positive correlation with tumor proximity")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_correlation_T_breast.pdf",height=7,width=13)
patchwork::wrap_plots(p1,p2,ncol=2,guides="collect")
dev.off()

# pdf("~/desktop/boundary_vs_GS.pdf",onefile = T,height=12,width=8)
# for(i in 1:28){
#   gs = tmp[,1:28][,i]
#   lm = lm(tmp[,i] ~ tmp$V1)
#   p9=ggplot(tmp, aes(x = V1, y = gs)) +
#     geom_point() +
#     geom_abline(slope = coef(lm)[2], intercept = coef(lm)[1])+
#     geom_smooth(span=0.1)+
#     xlab("")+ylab(names(tmp[,1:28])[i])+
#     annotate("text",x=0.7,y=max(tmp[,1:28][,i]),label=paste("Correlation:",round(cor.test(tmp$V1,tmp[,1:28][,i])$estimate,2)))+
#     annotate("text",x=0.7,y=max(tmp[,1:28][,i]-(max(tmp[,1:28][,i])-min(tmp[,1:28][,i]))/12),label=paste("p-value:",format.pval(cor.test(tmp$V1,tmp[,1:28][,i])$p.value)))+
#     theme_classic()
#   print(wrap_plots(p9,p1,p6,p4,p7,p8,ncol=1)+ plot_layout(heights = c(3,1,1,1,1,1)))
# }
# dev.off()


############################ Main Figure 6 ############################
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_Signature_score_G1.pdf",height=10, width=10,onefile = T)
pdf("~/desktop/HALLMARK_Spatial_Signature_score_G1.pdf",height=10, width=10,onefile = T)
for (i in 41:90){ 
  # colnames(brast_ST@meta.data)[which(colnames(brast_ST@meta.data)==colnames(brast_ST@meta.data)[i])] =
  #   substr(colnames(brast_ST@meta.data)[i],1,nchar(colnames(brast_ST@meta.data)[i])-3)
  print(SpatialFeaturePlot(brast_ST, features = colnames(brast_ST@meta.data)[i])+
          theme(text = element_text(size = 20, face = "bold"))+
          theme(legend.text=element_text(size=10, color="black"))
          # scale_fill_gradient(limits = c(range(brast_ST@meta.data[,i])[1], range(brast_ST@meta.data[,i])[2]),
          #                      values = round(as.vector(summary(brast_ST@meta.data[,i])[c(1,4,5,6)]),2),
          #                      low = "blue", high = "brown4")
        )
}
# round(as.vector(summary(brast_ST@meta.data[,i])[c(1,4,6)]),2)
dev.off()

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_TRM.pdf",
    height=8,width=8)
SpatialFeaturePlot(brast_ST,features=c("Tissue.Resident.Memory1"))
dev.off()

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_PossibleTRM.pdf",
    height=4,width=12)
SpatialFeaturePlot(brast_ST,features=c("CD8A","CD3D","TIGIT"))
dev.off()

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_MP_TRM.pdf",
    height=4,width=12)
SpatialFeaturePlot(brast_ST,features=c("L_Chemokine_Mediated_T_activation1","L_CellCycle_Immune_Response1","L_Lymphocyte_Activation1"))
dev.off()
# brast_ST <- AddMetaData(object = brast_ST, metadata = meta)
# brast_ST$Annotation = ifelse(brast_ST$Annotation=="","others",brast_ST$Annotation)
# Idents(brast_ST) = "Annotation"
# DotPlot(brast_ST, features = colnames(brast_ST@meta.data[, c(4:29)])) + RotatedAxis()


# ######################## scRNA-seq reference data ######################## 
# 
# metadata = as.data.frame(read.delim(paste0(path_data,"spatial/intestine-scRNA/intestine_meta.tsv"),row.names = 1))
# exp = Read10X_h5(paste0(path_data,"spatial/intestine-scRNA/intestine_expression.h5"), use.names = TRUE, unique.features = TRUE)
# exp = exp[,which(colnames(exp)%in%rownames(metadata))]
# # write.csv(as.matrix(exp),paste0(path_data,"spatial/intestine-scRNA/intestine_count.csv"))
# # write.csv(metadata,paste0(path_data,"spatial/intestine-scRNA/intestine_metadata.csv"))
# srt = CreateSeuratObject(exp)
# srt = AddMetaData(srt,metadata=metadata)
# srt <- SCTransform(srt, ncells = 3000, verbose = FALSE)
# srt <- RunPCA(srt, verbose = FALSE)
# srt <- RunUMAP(srt, dims = 1:30)
# saveRDS(srt,paste0(path_data,"spatial/intestine-scRNA/intestine.rds"))
# 
# metadata = as.data.frame(read.delim(paste0(path_data,"spatial/ovarian-scRNA/ovarian_meta.tsv"),row.names = 1))
# exp = Read10X_h5(paste0(path_data,"spatial/ovarian-scRNA/ovarian_expression.h5"), use.names = TRUE, unique.features = TRUE)
# exp = exp[,which(colnames(exp)%in%rownames(metadata))]
# # write.csv(as.matrix(exp),paste0(path_data,"spatial/ovarian-scRNA/ovarian_count.csv"))
# # write.csv(metadata,paste0(path_data,"spatial/ovarian-scRNA/ovarian_metadata.csv"))
# srt = CreateSeuratObject(exp)
# srt = AddMetaData(srt,metadata=metadata)
# srt <- SCTransform(srt, ncells = 3000, verbose = FALSE)
# srt <- RunPCA(srt, verbose = FALSE)
# srt <- RunUMAP(srt, dims = 1:30)
# saveRDS(srt,paste0(path_data,"spatial/ovarian-scRNA/ovarian.rds"))
# 
# metadata = as.data.frame(read.delim(paste0(path_data,"spatial/breast-scRNA/breast_meta.tsv"),row.names = 1))
# exp = Read10X_h5(paste0(path_data,"spatial/breast-scRNA/breast_expression.h5"), use.names = TRUE, unique.features = TRUE)
# exp = exp[,which(colnames(exp)%in%rownames(metadata))]
# # write.csv(as.matrix(exp),paste0(path_data,"spatial/breast-scRNA/breast_count.csv"))
# # write.csv(metadata,paste0(path_data,"spatial/breast-scRNA/breast_metadata.csv"))
# srt = CreateSeuratObject(exp)
# srt = AddMetaData(srt,metadata=metadata)
# srt <- SCTransform(srt, ncells = 3000, verbose = FALSE)
# srt <- RunPCA(srt, verbose = FALSE)
# srt <- RunUMAP(srt, dims = 1:30)
# saveRDS(srt,paste0(path_data,"spatial/breast-scRNA/breast.rds"))

# metadata = read.csv("/rsrch4/home/bcb/she4/hallmark/data/validation-data/GSE156728/cell.annotations.csv",row.names = 1)
# rownames(metadata) = metadata$cellID
# exp = read.csv("/rsrch4/home/bcb/she4/hallmark/data/validation-data/GSE156728/counts.csv")
# exp = exp[!duplicated(exp[,1]),]
# rownames(exp) = exp[,1]
# exp = exp[,-1]
# srt = CreateSeuratObject(as.matrix(exp))
# srt = AddMetaData(srt,metadata=metadata)
# srt <- SCTransform(srt, ncells = 3000, verbose = FALSE)
# srt <- RunPCA(srt, verbose = FALSE)
# srt <- RunUMAP(srt, dims = 1:30)
# saveRDS(srt,paste0(path_data,"spatial/T_reference.rds"))

# srt = readRDS(paste0(path_data,"spatial/breast-scRNA/breast.rds"))

######################## Cell Trek ######################## 

## Reference scRNA - remove the malignant cells
srt = readRDS(paste0(path_data,"spatial/breast-scRNA/breast.rds"))
# srt = subset(srt,subset = Celltype_minor%in%c("B-cells_Memory B-cells","B-cells_Naive B-cells",
#                                               "Myeloid_Cycling","Myeloid_DCs","Myeloid_Macrophage",
#                                               "Myeloid_Monocyte","T-cells_CD4+ T-cells","T-cells_CD8 T-cells",
#                                               "T-cells_Cycling","T-cells_NK cells","T-cells_NKT"))
#srt_subset = subset(srt,cells=colnames(srt)[which(srt$Celltype_major!="Cancer Epithelial")])
Idents(srt) = "Celltype_minor"
# DimPlot(srt_subset,group.by = "Celltype_minor",label=T)
srt_subset = subset(srt,downsample = 300)

## ST data to annotate - remove the malignant spots
# malignancy = read.delim(paste0(path_data,"spatial/breast/BC_copykat_prediction.txt"))
# rownames(malignancy) = malignancy$cell.names
# brast_ST = Load10X_Spatial(paste0(path_data,"spatial/breast/"))
# brast_ST = AddMetaData(brast_ST,metadata = malignancy)
# brast_ST_subset = subset(brast_ST,cells=colnames(brast_ST)[which(brast_ST$copykat.pred=="diploid")])
brast_ST_subset=brast_ST
traint <- CellTrek::traint(st_data=brast_ST_subset, sc_data=srt_subset, 
                                 sc_assay='RNA', cell_names='Celltype_major')
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=srt_subset, 
                                     sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)
saveRDS(celltrek,paste0(path_result,"spatial/celltrek_breast.rds"))
celltrek = readRDS(paste0(path_result,"spatial/celltrek_breast.rds"))
celltrek$cell_type <- factor(celltrek$celltrek$Celltype_major, levels=sort(unique(celltrek$celltrek$Celltype_major)))
celltrek1 = celltrek$celltrek
CellTrek::celltrek_vis(celltrek1@meta.data %>% dplyr::select(coord_x, coord_y, Celltype_major:id_new),
                       celltrek1@images$slice1@image,celltrek1@images$slice1@scale.factors$lowres)

# gs_list = c(gs_list,H_list,KEGG_list)
for (j in 1:length(gs_list[88:115])){
  marker_list <- gs_list[88:115][[j]][gs_list[88:115][[j]] %in% rownames(celltrek1@assays$RNA@data)]
  print (j); length(marker_list)
  name_new = names(gs_list[88:115])[j]
  celltrek1 <- AddModuleScore(celltrek1, features = list(marker_list), nbin = 12, name = name_new)
  names(celltrek1@meta.data)[which(names(celltrek1@meta.data)==paste0(name_new,"1"))] = name_new
}
celltrek1$Celltype_minor = factor(celltrek1$Celltype_minor,
                                  levels = c("Cancer Epithelial_Basal_SC","Cancer Epithelial_Her2_SC","Cancer Epithelial_LumA_SC","Cancer Epithelial_LumB_SC",
                                             "Cancer Epithelial_Cycling","Endothelial_Lymph_ECs","Endothelial_Endothelial ACKR1+","Endothelial_Endothelial RGS5+",
                                             "Endothelial_Endothelial CXCL12+","CAFs_Inflammatory-CAFs","CAFs_Myofibroblast-like CAFs",
                                             "SMCs_Immature SMCs","SMCs_Differeniated SMCs","SMCs_Cycling",
                                             "Myeloid_Cycling","Myeloid_Monocyte","Myeloid_Macrophage",
                                             "T-cells_CD8+ T-cells","T-cells_NK cells" ,"T-cells_CD4+ T-cells","T-cells_NKT","T-cells_Cycling",
                                             "Plasmablasts_Plasmablasts","B-cells_Naive B-cells","B-cells_Memory B-cells"))
Idents(celltrek1) = "Celltype_minor"
celltrek1 = subset(celltrek1,cells=colnames(celltrek1)[-which(is.na(celltrek1$Celltype_minor))])
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_clustered_dotplot.pdf",
#     height=9, width=11)
pdf("~/desktop/chloe_KEGG.pdf",height=7, width=10)
Clustered_DotPlot(celltrek1,features = names(gs_list),k=3, # [88:115]
                  colors_use_exp = c("white","wheat1","yellow","violetred","purple4"),
                  colors_use_idents = c(colorRampPalette(c("salmon", "brown"))(5),
                                        colorRampPalette(c("aquamarine", "darkgreen"))(4),
                                        "sienna2","tan3",
                                        colorRampPalette(c("azure2", "azure4"))(3),
                                        colorRampPalette(c("cadetblue1", "blue"))(3),
                                        colorRampPalette(c("plum2", "mediumpurple2"))(5),
                                        colorRampPalette(c("khaki1", "khaki4"))(3)))
dev.off()
print(a + theme(panel.grid.major = element_blank()))

celltrek1$Celltype_minor = factor(celltrek1$Celltype_minor,
                                 levels = c("Cancer Epithelial_Basal_SC","Cancer Epithelial_Her2_SC","Cancer Epithelial_LumA_SC","Cancer Epithelial_LumB_SC",
                                 "Cancer Epithelial_Cycling","Endothelial_Lymph_ECs","Endothelial_Endothelial ACKR1+","Endothelial_Endothelial RGS5+",
                                 "Endothelial_Endothelial CXCL12+","CAFs_Inflammatory-CAFs","CAFs_Myofibroblast-like CAFs",
                                 "SMCs_Immature SMCs","SMCs_Differeniated SMCs","SMCs_Cycling",
                                 "Myeloid_Cycling","Myeloid_Monocyte","Myeloid_Macrophage",
                                 "T-cells_CD8+ T-cells","T-cells_NK cells" ,"T-cells_CD4+ T-cells","T-cells_NKT","T-cells_Cycling",
                                 "Plasmablasts_Plasmablasts","B-cells_Naive B-cells","B-cells_Memory B-cells"))
Idents(celltrek1) = "Celltype_minor"
celltrek1 = subset(celltrek1,cells=colnames(celltrek1)[-which(is.na(celltrek1$Celltype_minor))])
Clustered_DotPlot(celltrek1,features = names(gs_list))
pdf("~/desktop/Spatial_dotplot_KEGG.pdf",height=10, width=15)
DotPlot(celltrek1,features = names(gs_list)[79:107])+coord_flip()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(text = element_text(face="bold"))
dev.off()
celltrek1$Tumor_vs_Immune = ifelse(celltrek1$Celltype_minor %in% c("Cancer Epithelial_Basal_SC","Cancer Epithelial_Her2_SC","Cancer Epithelial_LumA_SC","Cancer Epithelial_LumB_SC",
                                            "Cancer Epithelial_Cycling"),"Tumor",
                                   ifelse(celltrek1$Celltype_minor %in% c("Myeloid_Cycling","Myeloid_Monocyte","Myeloid_Macrophage",
                                                                         "T-cells_CD8+ T-cells","T-cells_NK cells" ,"T-cells_CD4+ T-cells","T-cells_NKT","T-cells_Cycling",
                                                                         "Plasmablasts_Plasmablasts","B-cells_Naive B-cells","B-cells_Memory B-cells"),"Immune",""))
celltrek2 = subset(celltrek1,cells=colnames(celltrek1)[-which(celltrek1$Tumor_vs_Immune=="")])
Idents(celltrek2) = "Tumor_vs_Immune"
pdf("~/desktop/Spatial_dotplot_tumor-immune_MP.pdf",height=10, width=7)
DotPlot(celltrek2,features = names(gs_list)[1:28])+coord_flip()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ 
  theme(text = element_text(face="bold"))
dev.off()

for (i in c("TCR_Anchoring","Interferon_induced_Antiviral_Defense",
            "IFN.gamma.induced_cytotoxicity","Chemokine_Mediated_immunity")){
  assign(paste0("p_",i),SpatialFeaturePlot(celltrek1, features = i,pt.size.factor = 1.2)+
          theme(text = element_text(size = 20, face = "bold"))+
          theme(legend.text=element_text(size=7, color="black")))
}
pdf("~/desktop/Spatial_Signature_score_4.pdf",height=7, width=18)
ggarrange(p_TCR_Anchoring,p_IFN.gamma.induced_cytotoxicity,p_Interferon_induced_Antiviral_Defense,
          nrow=1,ncol=3)
dev.off()

pdf("~/desktop/Spatial_Signature_score_4.pdf",height=7, width=7)
for (i in names(gs_list)){
  print(SpatialFeaturePlot(celltrek1, features = i,pt.size.factor = 1.2)+
    theme(text = element_text(size = 20, face = "bold"))+
    theme(legend.text=element_text(size=7, color="black")))
  print(i)
}
dev.off()

######################## Deconvolution ######################## 

# srt = readRDS(paste0(path_data,"spatial/T_reference.rds"))
# brast_ST = Load10X_Spatial(paste0(path_data,"spatial/breast/"))
# 
# plot1 <- VlnPlot(brast_ST, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(brast_ST, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)
# 
# brast_ST <- SCTransform(brast_ST, assay = "Spatial", verbose = FALSE)
# SpatialFeaturePlot(brast_ST, features = "nCount_Spatial") + theme(legend.position = "right")
# 
# brast_ST <- RunPCA(brast_ST, assay = "SCT", verbose = FALSE)
# brast_ST <- FindNeighbors(brast_ST, reduction = "pca", dims = 1:30)
# brast_ST <- FindClusters(brast_ST, verbose = FALSE)
# brast_ST <- RunUMAP(brast_ST, reduction = "pca", dims = 1:30)
# p1 <- DimPlot(brast_ST, reduction = "umap", label = TRUE)
# p2 <- SpatialDimPlot(brast_ST, label = TRUE, label.size = 3)
# wrap_plots(p1, p2)
# # DimPlot(srt, group.by = "Celltype..minor.lineage.", label = TRUE,raster=FALSE)
# 
# ##### Deconvolution by seurat default
# anchors <- FindTransferAnchors(reference = srt, query = brast_ST, normalization.method = "SCT")
# predictions.assay <- TransferData(anchorset = anchors, refdata = srt$cell.types, prediction.assay = TRUE,
#                                   weight.reduction = brast_ST[["pca"]], dims = 1:30)
# brast_ST[["predictions"]] <- predictions.assay
# DefaultAssay(brast_ST) <- "predictions"
# pdf(paste0(path_data,"spatial/breast/breast_spatial_sc_deconvoluted_T.pdf"),onefile = T)
# for (i in unique(srt$cell.types)){
#  print(SpatialFeaturePlot(brast_ST,
#                           features = i,
#                           pt.size.factor = 1.6)) 
# }
# dev.off()
# saveRDS(sbrast_ST,paste0(path_data,"breast/breast_spatial_deconvoluted.rds"))


###### Deconvolution by RCTD

# Idents(srt) <- "celltype_minor"
# # extract information to pass to the RCTD Reference function
# counts <- srt[["RNA"]]@counts
# cluster <- as.factor(srt$celltype_minor)
# names(cluster) <- colnames(srt)
# nUMI <- srt$nCount_RNA
# names(nUMI) <- colnames(srt)
# reference <- Reference(counts, cluster, nUMI)
# set up query with the RCTD function SpatialRNA
# counts <- brast_ST[["Spatial"]]@counts
# coords <- GetTissueCoordinates(brast_ST)
# colnames(coords) <- c("x", "y")
# coords[is.na(colnames(coords))] <- NULL
# query <- SpatialRNA(coords, counts, colSums(counts))
# 
# RCTD <- create.RCTD(query, reference, max_cores = 8)
# RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# brast_ST <- AddMetaData(brast_ST, metadata = RCTD@results$results_df)
# saveRDS(brast_ST,"/rsrch4/home/bcb/she4/hallmark/data/Wu_etal_2021_BRCA_scRNASeq/breast_spatial_deconvoluted.rds")
# 
# brast_ST = readRDS("/Volumes/she4/hallmark/data/Wu_etal_2021_BRCA_scRNASeq/breast_spatial_deconvoluted.rds")
# DefaultAssay(brast_ST) = "Spatial"
# SpatialDimPlot(brast_ST, group.by = "first_type")

# p1 <- SpatialDimPlot(brast_ST, group.by = "first_type")
# p2 <- SpatialDimPlot(brast_ST, group.by = "second_type")
# p1 | p2


# genes = read.delim("/rsrch4/home/bcb/she4/hallmark/data/spatial/intestine-scRNA/matrix.genes.tsv",header=F)
# cells = read.delim("/rsrch4/home/bcb/she4/hallmark/data/spatial/intestine-scRNA/matrix.barcodes.tsv",header=F)
# rownames(mat) = genes$V1
# colnames(mat) = cells$V1
# 
# meta = read.delim("/rsrch4/home/bcb/she4/hallmark/data/spatial/intestine-scRNA/metatable.tsv",header = F,row.names = 1)
# meta = meta[-2,]
# names(meta) = meta[1,]
# meta = meta[-1,]
