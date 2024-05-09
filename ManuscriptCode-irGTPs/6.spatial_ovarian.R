########## (0_ovarian) Create inputs for CytoSpace ##########
files = list.files("/Volumes/she4/hallmark/data/spatial/ovarian-scRNA/",pattern = "matrix.mtx",full.names = T)
barcode_files = list.files("/Volumes/she4/hallmark/data/spatial/ovarian-scRNA/",pattern = "barcode",full.names = T)
feature_files = list.files("/Volumes/she4/hallmark/data/spatial/ovarian-scRNA/",pattern = "genes",full.names = T)
mat_all = readMM(files[1])
colnames(mat_all) = read.table(barcode_files[1])$V1
rownames(mat_all) = read.table(feature_files[1])$V2
for(i in 2:length(files)){
  mat = readMM(files[i])
  colnames(mat) = read.table(barcode_files[i])$V1
  rownames(mat) = read.table(feature_files[i])$V2
  mat_all = cbind(mat_all,mat)
  print(i)
}
ov_scRNA_celllabels = read.csv("/Volumes/she4/hallmark/data/spatial/ovarian-scRNA/OV_GSE154600_CellMetainfo_table.tsv",sep="")
ov_scRNA_celllabels$Cell = unlist(strsplit(ov_scRNA_celllabels$Cell,"@"))[seq(2,nrow(ov_scRNA_celllabels)*2,2)]
ov_scRNA_celllabels = ov_scRNA_celllabels[,c("Cell","X.major.lineage.")]
names(ov_scRNA_celllabels) = c("SpotID","CellType")
common_cells = intersect(ov_scRNA_celllabels$SpotID,colnames(mat_all))
ov_scRNA_celllabels = ov_scRNA_celllabels[!duplicated(ov_scRNA_celllabels$SpotID),]
ov_scRNA_celllabels = ov_scRNA_celllabels[which(ov_scRNA_celllabels$SpotID%in%common_cells),]
mat_all = mat_all[,common_cells]
# mat_all = as.data.frame(cbind(rownames(mat_all),as.matrix(mat_all)))
# names(mat_all)[1] = "Gene"
scrna_count <- as.data.frame(as.matrix(mat_all))
scrna_count <- cbind(rownames(mat_all), scrna_count)
colnames(scrna_count)[1] <- 'GENES'
write.table(ov_scRNA_celllabels,"/Volumes/she4/hallmark/data/spatial/CytoSPACE_ovarian/ov_scRNA_celllabels.txt",
            row.names = F,quote = F,sep="\t")
write.table(scrna_count,"/Volumes/she4/hallmark/data/spatial/CytoSPACE_ovarian/ov_scRNA_GEP.txt",
          row.names = T,quote = F,sep="\t")

ovarian = Load10X_Spatial(paste0(path_data,"spatial/ovarian/"))
ov_STdata_coordinates = cbind(colnames(ovarian),ovarian@images$slice1@coordinates[,c("row","col")])
names(ov_STdata_coordinates)[1] = "SpotID"
write.table(ov_STdata_coordinates,"/Volumes/she4/hallmark/data/spatial/CytoSPACE_ovarian/ov_STdata_coordinates.txt",
            row.names = F,quote = F,sep="\t")
ov_STdata_GEP = as.data.frame(cbind(rownames(ovarian),as.matrix(ovarian@assays$Spatial$counts)))
names(ov_STdata_GEP)[1] = "Gene"
write.table(ov_STdata_GEP,"/Volumes/she4/hallmark/data/spatial/CytoSPACE_ovarian/ov_STdata_GEP.txt",
            row.names = F,quote = F,sep="\t")


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
effector=c('NKG7','PRF1','GZMA','GZMB','GZMH','GNLY')
exh = read.delim(paste0(path_data,"exh.txt"),header = F)
exh = exh$V1
naive=c('CCR7','LEF1','SELL','TCF7') #from CXCR6/CXCL16 paper
tissue_memory.full=c('ITGAE','ITGA1','ZNF683','CD69','IL2','IL10','CA10','CXCR6','CXCL13','DUSP6','KCNK5','RGS1','CRTAM','PDCD1')
TEM = c("CX3CR1","GZMB","GNLY")
CD4 = c("CD4")
CD8 = c("CD8A","CD69")
Treg = c("FOXP3","IL2RA")
M1 = c("CD80","CD86","CD68","IL1R1","TLR2","TLR4","SOCS3","TNF")
M2 = c("CD163","MRC1","CD200R1","TGM2","IL1R2","TGFB1","IL10")
genesetlist=list(exh,exhaustion,effector,naive,tissue_memory.full,TEM,CD4,CD8,Treg,M1,M2)
genesetlist_names=c("Exhaustion",'exhaustion_short','Cytotoxicity',
                    'naive','Tissue.Resident.Memory',"Tissue.Effector.Memory","CD4","CD8","Treg","M1",
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

####################### Supplementary to Figure 6 #######################
brast_ST = Load10X_Spatial(paste0(path_data,"spatial/ovarian/"))
for (j in 1:length(gs_list)){
  marker_list <- gs_list[[j]][gs_list[[j]] %in% rownames(brast_ST@assays$Spatial@counts)]
  print (j); length(marker_list)
  name_new = names(gs_list)[j]
  brast_ST <- AddModuleScore(brast_ST, features = list(marker_list), nbin = 12, name = name_new)
}
TLS = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/10X_pathologist/ovarian/TLS.csv",
               row.names = 1)
TLS$TLS = ifelse(TLS$TLS=="TLS","TLS","Not TLS")
brast_ST = AddMetaData(brast_ST,metadata = TLS)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_TLS_ovarian.pdf",
    height=4,width=4)
VlnPlot(brast_ST,features = "L_TCR_Anchoring1",group.by = "TLS",pt.size=0,
        cols = c("hotpink","cornflowerblue"))+
  ggpubr::stat_compare_means(label="p.signif")+ylim(-2.5,3)
dev.off()
malignancy = read.delim(paste0(path_data,"spatial/ovarian/OC_copykat_prediction.txt"))
rownames(malignancy) = malignancy$cell.names
brast_ST = AddMetaData(brast_ST,metadata = malignancy)


######## Incorporate CytoSpace celltype deconvolution ########
# cell_type_assignments_by_spot = read.csv("/Volumes/she4/hallmark/results/spatial/cytospace_results/cell_type_assignments_by_spot.csv",row.names = 1)
# cell_type_assignments_by_spot$T_percent = (cell_type_assignments_by_spot$T.cells.CD8+cell_type_assignments_by_spot$T.cells.CD4)/cell_type_assignments_by_spot$Total.cells
# cell_type_assignments_by_spot$B_percent = (cell_type_assignments_by_spot$B.cells)/cell_type_assignments_by_spot$Total.cells
# cell_type_assignments_by_spot$Myeloid_percent = (cell_type_assignments_by_spot$Monocytes.and.Macrophages)/cell_type_assignments_by_spot$Total.cells
# cell_type_assignments_by_spot$barcode = unlist(strsplit(rownames(cell_type_assignments_by_spot),"__"))[seq(1,nrow(cell_type_assignments_by_spot)*2,2)]
# cell_type_assignments_by_spot$barcode = gsub("\\.","-",cell_type_assignments_by_spot$barcode)
# rownames(cell_type_assignments_by_spot) = cell_type_assignments_by_spot$barcode
# brast_ST = AddMetaData(brast_ST,metadata = cell_type_assignments_by_spot)

########## Tumor Boundaries and gene set relationship ##########
brast_ST = UpdateSeuratForSemla(object=brast_ST,verbose = TRUE)
brast_ST <- SetIdent(brast_ST, value = "cell_names")
for(s in 1:length(unique(brast_ST$cell.names))){ # zero neighbors: 1715
  # if(s == 1715)next
  brast_ST <- RegionNeighbors(brast_ST, column_name = "cell.names",
                              column_labels = unname(unique(brast_ST$cell.names)[s]), 
                              verbose = TRUE)
  print(s)
}
spot = names(brast_ST@meta.data)[140:ncol(brast_ST@meta.data)]
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

brast_ST_T = subset(brast_ST,subset=Myeloid_percent>0.8&copykat.pred=="diploid")
tmp = brast_ST_T@meta.data[,c(4:92,which(names(brast_ST_T@meta.data)=="V1"))]

# tmp = brast_ST_T@meta.data[-which(brast_ST_T@meta.data$V1==1|brast_ST_T@meta.data$copykat.pred=="aneuploid"),c(4:90,which(names(brast_ST_T@meta.data)=="V1"))]

signatures = c("HALLMARK_ALLOGRAFT_REJECTION1","HALLMARK_INFLAMMATORY_RESPONSE1",
               "HALLMARK_INTERFERON_GAMMA_RESPONSE1",
               "naive1","Treg1","Tissue.Resident.Memory1",
               "Cytotoxicity1","Tissue.Effector.Memory1","Exhaustion1")
signatures = c("HALLMARK_ALLOGRAFT_REJECTION1","HALLMARK_INFLAMMATORY_RESPONSE1",
               "HALLMARK_INTERFERON_GAMMA_RESPONSE1","HALLMARK_INTERFERON_ALPHA_RESPONSE1",
               "M11","M21")
temp = tmp[,c("V1",signatures)]
temp$V1 = round(temp$V1,2)
temp$V1 = factor(temp$V1,sort(unique(temp$V1)))
temp = aggregate(temp[,signatures],by=list(temp$V1),FUN=mean)
temp_long = melt(setDT(temp),id.var="Group.1")
temp_long$value = round(temp_long$value,2)
temp_long$variable = factor(temp_long$variable,signatures)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_Infiltration_vs_Distance_Myeloid.pdf",
    height=4,width=11)
ggplot(temp_long,aes(x=Group.1,y=value))+
  geom_histogram(stat = "identity",position = position_dodge())+
  facet_wrap(~variable,scale="free_y",ncol=3)+
  xlab("Neighboring Malignant %")+ylab("Z-norm Cell Type Score")+
  geom_hline(yintercept=0, linetype="dashed", color = "brown3", size=0.6)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

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

# T
p1=ggpairs(tmp, columns = c(5,14,15,29,31,34),
           lower = list(continuous = wrap(lowerfun)))+
  theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))
p2=ggpairs(tmp, columns = c(4,17,2,32,33,37),lower = list(continuous = wrap(lowerfun)))+
  theme(strip.text.x = element_text(size = 8),strip.text.y = element_text(size = 8))


############################ Main Figure 6 ############################
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Spatial_Signature_score_G1.pdf",height=10, width=10,onefile = T)
pdf("~/desktop/MP_Spatial_Signature_score_G1_ovarian.pdf",height=10, width=10,onefile = T)
for (i in 4:31){ 
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


######################## Cell Trek ######################## 

srt = readRDS(paste0(path_data,"spatial/ovarian-scRNA/ovarian.rds"))
Idents(srt) = "Celltype..major.lineage."
srt_subset = subset(srt,downsample = 300)
brast_ST_subset
brast_ST_subset=brast_ST

traint <- CellTrek::traint(st_data=brast_ST_subset, sc_data=srt_subset, st_assay = "Spatial",
                           sc_assay='RNA', cell_names='Celltype..major.lineage.')
celltrek <- CellTrek::celltrek(st_sc_int=traint, int_assay='traint', sc_data=srt_subset, 
                               sc_assay = 'RNA', 
                               reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                               dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)
saveRDS(celltrek,paste0(path_result,"spatial/celltrek_ovarian.rds"))
celltrek = readRDS(paste0(path_result,"spatial/celltrek_ovarian.rds"))
celltrek$cell_type <- factor(celltrek$celltrek$Celltype..major.lineage., levels=sort(unique(celltrek$celltrek$Celltype..major.lineage.)))
celltrek1 = celltrek$celltrek
CellTrek::celltrek_vis(celltrek1@meta.data %>% dplyr::select(coord_x, coord_y, Celltype..major.lineage.:id_new),
                       celltrek1@images$slice1@image,celltrek1@images$slice1@scale.factors$lowres)

# gs_list = c(gs_list,H_list,KEGG_list)
for (j in 1:length(gs_list[1:118])){
  marker_list <- gs_list[1:118][[j]][gs_list[1:118][[j]] %in% rownames(celltrek1@assays$RNA@data)]
  print (j); length(marker_list)
  name_new = names(gs_list[1:118])[j]
  celltrek1 <- AddModuleScore(celltrek1, features = list(marker_list), nbin = 12, name = name_new)
  names(celltrek1@meta.data)[which(names(celltrek1@meta.data)==paste0(name_new,"1"))] = name_new
}
celltrek1$Celltype..major.lineage. = factor(celltrek1$Celltype..major.lineage.,
                                  levels = c("B","Plasma","CD8Tex","CD8T","Tprolif","DC","Mono/Macro","Mast",
                                             "Malignant","Epithelial","Fibroblasts","Endothelial","Myofibroblasts"))
Idents(celltrek1) = "Celltype..major.lineage."
celltrek1 = subset(celltrek1,cells=colnames(celltrek1)[-which(is.na(celltrek1$Celltype..major.lineage.))])
saveRDS(celltrek1,paste0(path_result,"spatial/celltrek_ovarian_cleaned.rds"))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/spatial_clustered_dotplot_ovarian_Hallmark.pdf",
    height=7, width=8)
Clustered_DotPlot(celltrek1,features = names(gs_list)[40:89],k=3, # [88:115]
                  colors_use_exp = c("white","wheat1","yellow","violetred","purple4"),
                  # colors_use_idents = c(colorRampPalette(c("salmon", "brown"))(5),
                  #                       colorRampPalette(c("aquamarine", "darkgreen"))(4),
                  #                       "sienna2","tan3",
                  #                       colorRampPalette(c("azure2", "azure4"))(3),
                  #                       colorRampPalette(c("cadetblue1", "blue"))(3),
                  #                       colorRampPalette(c("plum2", "mediumpurple2"))(5),
                  #                       colorRampPalette(c("khaki1", "khaki4"))(3)),
                  colors_use_idents = c("khaki1","khaki4","plum2", "mediumpurple2","cadetblue1",
                                                "azure4","salmon","brown","sienna2","tan3"),
                                                cluster_ident = F
                  )
dev.off()

# print(a + theme(panel.grid.major = element_blank()))
# 
# celltrek1$Celltype_minor = factor(celltrek1$Celltype_minor,
#                                   levels = c("Cancer Epithelial_Basal_SC","Cancer Epithelial_Her2_SC","Cancer Epithelial_LumA_SC","Cancer Epithelial_LumB_SC",
#                                              "Cancer Epithelial_Cycling","Endothelial_Lymph_ECs","Endothelial_Endothelial ACKR1+","Endothelial_Endothelial RGS5+",
#                                              "Endothelial_Endothelial CXCL12+","CAFs_Inflammatory-CAFs","CAFs_Myofibroblast-like CAFs",
#                                              "SMCs_Immature SMCs","SMCs_Differeniated SMCs","SMCs_Cycling",
#                                              "Myeloid_Cycling","Myeloid_Monocyte","Myeloid_Macrophage",
#                                              "T-cells_CD8+ T-cells","T-cells_NK cells" ,"T-cells_CD4+ T-cells","T-cells_NKT","T-cells_Cycling",
#                                              "Plasmablasts_Plasmablasts","B-cells_Naive B-cells","B-cells_Memory B-cells"))
# Idents(celltrek1) = "Celltype_minor"
# celltrek1 = subset(celltrek1,cells=colnames(celltrek1)[-which(is.na(celltrek1$Celltype_minor))])
# Clustered_DotPlot(celltrek1,features = names(gs_list))
# pdf("~/desktop/Spatial_dotplot_KEGG.pdf",height=10, width=15)
# DotPlot(celltrek1,features = names(gs_list)[79:107])+coord_flip()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
#   theme(text = element_text(face="bold"))
# dev.off()
# celltrek1$Tumor_vs_Immune = ifelse(celltrek1$Celltype_minor %in% c("Cancer Epithelial_Basal_SC","Cancer Epithelial_Her2_SC","Cancer Epithelial_LumA_SC","Cancer Epithelial_LumB_SC",
#                                                                    "Cancer Epithelial_Cycling"),"Tumor",
#                                    ifelse(celltrek1$Celltype_minor %in% c("Myeloid_Cycling","Myeloid_Monocyte","Myeloid_Macrophage",
#                                                                           "T-cells_CD8+ T-cells","T-cells_NK cells" ,"T-cells_CD4+ T-cells","T-cells_NKT","T-cells_Cycling",
#                                                                           "Plasmablasts_Plasmablasts","B-cells_Naive B-cells","B-cells_Memory B-cells"),"Immune",""))
# celltrek2 = subset(celltrek1,cells=colnames(celltrek1)[-which(celltrek1$Tumor_vs_Immune=="")])
# Idents(celltrek2) = "Tumor_vs_Immune"
# pdf("~/desktop/Spatial_dotplot_tumor-immune_MP.pdf",height=10, width=7)
# DotPlot(celltrek2,features = names(gs_list)[1:28])+coord_flip()+ 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("")+ 
#   theme(text = element_text(face="bold"))
# dev.off()
# 
# for (i in c("TCR_Anchoring","Interferon_induced_Antiviral_Defense",
#             "IFN.gamma.induced_cytotoxicity","Chemokine_Mediated_immunity")){
#   assign(paste0("p_",i),SpatialFeaturePlot(celltrek1, features = i,pt.size.factor = 1.2)+
#            theme(text = element_text(size = 20, face = "bold"))+
#            theme(legend.text=element_text(size=7, color="black")))
# }
# pdf("~/desktop/Spatial_Signature_score_4.pdf",height=7, width=18)
# ggarrange(p_TCR_Anchoring,p_IFN.gamma.induced_cytotoxicity,p_Interferon_induced_Antiviral_Defense,
#           nrow=1,ncol=3)
# dev.off()
# 
# pdf("~/desktop/Spatial_Signature_score_4.pdf",height=7, width=7)
# for (i in names(gs_list)){
#   print(SpatialFeaturePlot(celltrek1, features = i,pt.size.factor = 1.2)+
#           theme(text = element_text(size = 20, face = "bold"))+
#           theme(legend.text=element_text(size=7, color="black")))
#   print(i)
# }
# dev.off()
