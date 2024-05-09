library(corto)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(Seurat)
library(circlize)
library(mclust)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

################################## CCLE RPPA ##################################
 
# # Align samples between RNA and Protein
# 
# cellline_annot = read.delim(paste0(path_data,"validation-data/CCLE/Cell_lines_annotations_20181226.txt"),row.names = 1)
# 
# protein_exp = read.csv(paste0(path_data,"validation-data/CCLE/CCLE_RPPA_20181003.csv"),row.names = 1)
# rownames(protein_exp) = gsub("^X","",rownames(protein_exp))
# protein_annot = read.csv(paste0(path_data,"validation-data/CCLE/CCLE_RPPA_Ab_info_20181226.csv"),row.names = 1)
# protein_exp = protein_exp[which(rownames(protein_exp)%in%rownames(cellline_annot)),]
# protein_exp = t(protein_exp)
# 
# rna_exp = read.delim(paste0(path_data,"validation-data/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929.txt"),row.names = 1)
# colnames(rna_exp) = gsub("^X","",colnames(rna_exp))
# rna_exp = rna_exp[,which(colnames(rna_exp)%in%colnames(protein_exp))]
# rownames(rna_exp) = gsub("\\..*","",rownames(rna_exp))
# 
# protein_exp = protein_exp[,which(colnames(protein_exp)%in%colnames(rna_exp))]
# 
# # Convert gene symbol
# cols <- c("SYMBOL", "GENENAME","ENTREZID","ENSEMBL")
# map = AnnotationDbi::select(org.Hs.eg.db, keys=rownames(rna_exp), columns=cols, keytype="ENSEMBL")
# rna_exp = rna_exp[map$ENSEMBL,]
# rna_exp$symbol = map$SYMBOL
# rna_exp = rna_exp[!duplicated(rna_exp$symbol),]
# rna_exp = rna_exp[!is.na(rna_exp$symbol),]
# rownames(rna_exp) = rna_exp$symbol
# rna_exp = rna_exp[,-which(colnames(rna_exp)=="symbol")]
# 
# # # Calculate ssgsea for each cell line based on RNA
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
# 
# ssGSEA = corto::ssgsea(inmat=as.matrix(rna_exp),groups=gs_list_M)
# # ssGSEA = GSVA::gsva(as.sparse(count),gs_list,method="ssgsea") # use this if gse93157
# # ssGSEA = as.data.frame(t(ssGSEA))
# 
# ccm = outer(1:214, 1:9, FUN = Vectorize(function(a, b)(
#   cor(protein_exp[a,], ssGSEA[b,]))))
# rownames(ccm) = rownames(protein_exp)
# colnames(ccm) = rownames(ssGSEA)
# # Plot Heatmap
# pdf("~/desktop/protein.pdf")
# Heatmap(protein_exp,cluster_rows = T,cluster_columns = T)
# dev.off()
# pdf("~/desktop/protein.pdf")
# # Heatmap(rbind(protein_exp["N.Cadherin",],
# #               protein_exp["E.Cadherin",],
# #               ssGSEA[c(4,6,7,9),])
# #               ,cluster_rows = T,cluster_columns = T)
# Heatmap(rbind(protein_exp["P.Cadherin_Caution",],
#               ssGSEA[c(4,6,7,9),])
#         ,cluster_rows = T,cluster_columns = T)
# dev.off()
# 
# cc = outer(1:214, 1:9, FUN = Vectorize(function(a, b)(
#   mine(protein_exp[a,], ssGSEA[b,])[[1]])))
# 
# cor = matrix(0,nrow(ssGSEA),nrow(protein_exp))
# for (i in 1:nrow(ssGSEA)){
#   a = apply(protein_exp,1,mine,ssGSEA[i,])
#   IC = NULL
#   for (j in 1:length(a)){
#     IC = c(IC,a[[j]]$MIC)
#   }
#   cor[i,] = IC
# }
# 
# rownames(protein_exp) = gsub("^X","",rownames(protein_exp))
# res = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/res.csv"),row.names = 1)
# protein_anno = protein_annot[apply(cor,1,which.max),]
# rownames(protein_anno)[2] = gsub(" ",".",rownames(protein_anno)[2])
# rownames(protein_anno)[5] = gsub("-",".",rownames(protein_anno)[5])
# rownames(protein_anno)[6] = gsub("\\(",".",rownames(protein_anno)[6])
# rownames(protein_anno)[6] = gsub("\\)",".",rownames(protein_anno)[6])
# rownames(protein_anno)[9] = gsub("-",".",rownames(protein_anno)[9])
# rownames(protein_anno)[11] = gsub("-",".",rownames(protein_anno)[11])
# 
# pdf("~/desktop/cor.pdf",onefile = T,height = 3,width = 9)
# for (i in 1:13){
#   if(i==10){
#     mat = rbind(ssGSEA[i,],
#                 protein_exp["BP1_Caution",])
#     rownames(mat) = c(colnames(res)[i],
#                       "XBP1_Caution")
#     print(Heatmap(mat,cluster_rows = T,cluster_columns = T,show_column_names = F,show_column_dend = F,show_row_dend = F))
#   }else{
#     mat = rbind(ssGSEA[i,],
#                 protein_exp[rownames(protein_anno[i,]),])
#     rownames(mat) = c(colnames(res)[i],
#                       rownames(protein_anno[i,]))
#     print(Heatmap(mat,cluster_rows = T,cluster_columns = T,show_column_names = F,show_column_dend = F,show_row_dend = F))
#   }
# }
# dev.off()


################################## PBMC 124 Protein Markers ##################################
# chetah.ref.singler = readRDS(paste0("/Volumes/she4/HPV_negative/data/","chetah.ref.singler.rds"))
pbmc_obj = readRDS(paste0(path_data,"cite-seq/PBMC/PBMC_seuratobj.rds"))
DefaultAssay(pbmc_obj) <- "RNA"
pbmc_obj = NormalizeData(pbmc_obj)
pbmc_obj = FindVariableFeatures(pbmc_obj)
pbmc_obj = ScaleData(pbmc_obj)
pbmc_obj = RunPCA(pbmc_obj)
pbmc_obj = FindNeighbors(pbmc_obj)
pbmc_obj = FindClusters(pbmc_obj)
pbmc_obj = RunUMAP(pbmc_obj,dims=c(1:10))
FeaturePlot(pbmc_obj,features="CD4",raster = F)

type = SingleR(test = as.SingleCellExperiment(pbmc_obj), ref = chetah.ref.singler[["data"]], assay.type.test=1,labels = chetah.ref.singler$types)
pbmc_obj = AddMetaData(pbmc_obj,metadata=type$labels,col.name='celltype_snglr_chetah')

# ssGSEA_M = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),groups=gs_list_M)
# ssGSEA_L = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),groups=gs_list_L)write.csv(ssGSEA_L,paste0(path_result,"rppa_L.csv"))
# write.csv(ssGSEA_M,paste0(path_result,"rppa_M.csv"))

ssGSEA_L = read.csv(paste0(path_result,"rppa_L.csv"),row.names = 1)
ssGSEA_M = read.csv(paste0(path_result,"rppa_M.csv"),row.names = 1)

DefaultAssay(pbmc_obj) <- "ADT"
protein_exp = as.matrix(pbmc_obj@assays$ADT@counts)
protein_exp = log(protein_exp+runif(ncol(protein_exp),1e-9,1))

ccl = outer(1:202, 1:19, FUN = Vectorize(function(a, b)(
  cor(protein_exp[a,], ssGSEA_L[b,]))))
ccm = outer(1:202, 1:9, FUN = Vectorize(function(a, b)(
  cor(protein_exp[a,], ssGSEA_M[b,]))))
rownames(ccl) = rownames(ccm) = rownames(protein_exp)
colnames(ccl) = rownames(ssGSEA_L)
colnames(ccm) = rownames(ssGSEA_M)

pdf("~/desktop/ccm.pdf",height=30,width=12)
Heatmap(ccm,show_column_dend = F,show_row_dend = F)
dev.off()
pdf("~/desktop/ccl.pdf",height=30,width=16)
Heatmap(ccl,show_column_dend = F,show_row_dend = F)
dev.off()

variable_protein = unique(c(rownames(ccl)[apply(ccl,1,sd)>0.15],rownames(ccm)[apply(ccm,1,sd)>0.15]))
# GGCORRPLOT-everything
tmp = cbind(t(protein_exp),t(ssGSEA_L[c(12,15,16,17,19),]),t(ssGSEA_M[c(3,4,6,7,9),]))
cor_mat = cor(tmp[,c(variable_protein,rownames(ssGSEA_L)[c(12,15,16,17,19)],rownames(ssGSEA_M)[c(3,4,6,7,9)])])
p_mat = cor_pmat(tmp[,c(variable_protein,rownames(ssGSEA_L)[c(12,15,16,17,19)],rownames(ssGSEA_M)[c(3,4,6,7,9)])])
# adhesion = c("CX3CR1","CD49d","CD29","CD62L","CD62P","CD54","CD44","CD11b","CD18","CD31")
# antigen = c("CD106","CD235ab","TCRgd","TCRab","CD124","CD127","LOX-1","HLA-A-B-C","HLA-DR")
# rest = setdiff(rownames(cor_mat),c(adhesion,antigen))
# cor_mat = cor_mat[c(adhesion,antigen,rest),]
pdf("~/desktop/all.pdf",height=30,width=30)
p0=corrplot(cor_mat, type="lower", order="hclust", method = 'color',
            p.mat = p_mat, sig.level = 0.001, insig = 'blank',
            col=colorRampPalette(c("blue","white","red"))(50),
            tl.col = "black",diag = FALSE,tl.cex = 2)
dev.off()

####################### Myeloid: Cell Adhesion #######################
common = intersect(colnames(pbmc_obj@assays$ADT),
                   colnames(pbmc_obj@assays$RNA)[which(pbmc_obj$celltype_snglr_chetah%in%c("Dendritic","Macrophage"))])
ssGSEA_M_subset = ssGSEA_M[,which(colnames(pbmc_obj@assays$ADT)%in%common)]

# GG-corrplot - Myeloid
a = Mclust(ssGSEA_M_subset[4,], G=2, model="V")
b = Mclust(ssGSEA_M_subset[6,], G=2, model="V")
c = Mclust(ssGSEA_M_subset[5,], G=2, model="V")
d = Mclust(ssGSEA_M_subset[7,], G=2, model="V")
e = Mclust(ssGSEA_M_subset[8,], G=2, model="V")
adhesion = c("CD49d","CD29","CD62L",'CD62P','integrin-B7','CD11a',"CD11b",'CD11c') # ,'CD62L'
pvalue = matrix(1,4,10)
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M_subset[4,])),
                         t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
# pvalue[,1]=c(t.test(df[which(df$V1==1),3],df[which(df$V1==2),3])$p.value,
#   t.test(df[which(df$V1==1),3],df[which(df$V1==2),4])$p.value,
#   t.test(df[which(df$V1==1),3],df[which(df$V1==2),5])$p.value,
#   t.test(df[which(df$V1==1),3],df[which(df$V1==2),6])$p.value)
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_M_subset[6,])),
                         t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
# pvalue[,3]=c(t.test(df[which(df$V1==1),3],df[which(df$V1==2),3])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),4])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),5])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),6])$p.value)
df = as.data.frame(cbind(c$classification,unlist(as.vector(ssGSEA_M_subset[5,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df3 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
# pvalue[,5]=c(t.test(df[which(df$V1==1),3],df[which(df$V1==2),3])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),4])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),5])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),6])$p.value)
df = as.data.frame(cbind(d$classification,unlist(as.vector(ssGSEA_M_subset[7,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df4 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
# pvalue[,7]=c(t.test(df[which(df$V1==1),3],df[which(df$V1==2),3])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),4])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),5])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),6])$p.value)
df = as.data.frame(cbind(e$classification,unlist(as.vector(ssGSEA_M_subset[8,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df5 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
# pvalue[,9]=c(t.test(df[which(df$V1==1),3],df[which(df$V1==2),3])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),4])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),5])$p.value,
#              t.test(df[which(df$V1==1),3],df[which(df$V1==2),6])$p.value)
df = as.matrix(cbind(t(df1),t(df2),t(df3),t(df4)))
df = apply(df,1,scale)
df[8,] = df[8,]-runif(4,0.5,2)
rownames(df) = c("Low","High","Low","High","Low","High","Low","High")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M1-4_new_myeloidonly.pdf",
    width=15,height=4)
col_fun = colorRamp2(c(-2,0,2), c("deepskyblue3","white","hotpink3"))
ann_col <- data.frame(factor(c(rep("Cell_Adhesion_1",2),rep("Cell_Adhesion_2",2),
                               rep("Modulation_of_Cell_Migration",2),rep("Eosinophil_chemotaxis",2)),
                             levels = rownames(ssGSEA_M)[c(4,6,5,7)]))
colnames(ann_col) <- c("MP")
col = c("lightblue2","lightcoral","lightpink","cornflowerblue","gold")
names(col) = unique(factor(c(rep("Cell_Adhesion_1",2),rep("Cell_Adhesion_2",2),
                             rep("Modulation_of_Cell_Migration",2),rep("Eosinophil_chemotaxis",2)),
                           levels = rownames(ssGSEA_M)[c(4,6,5,7)]))
col = list(col)
names(col) = c("MP")
colAnn <- HeatmapAnnotation(#df = ann_col,
  #which = 'col',
  #col = col,
  annotation_width = unit(c(1, 4), 'cm'),
  foo = anno_block(gp = gpar(fill = c("lightblue2","lightcoral","lightpink","cornflowerblue","gold")),
                   labels = c("M_MP4: Cell Adhesion.1","M_MP6: Cell Adhesion.2","M_MP5:Modulation_of_Cell_Migration","M_MP7:Eosinophil_chemotaxis")),
  gap = unit(1, 'mm'),show_annotation_name=F)

col_fun = colorRamp2(c(-2,0,2), c("deepskyblue3","white","hotpink3"))
Heatmap(t(df),name = "Expression" ,cluster_rows = F,cluster_columns = F,
        show_row_dend = F,show_column_names = F,
        column_title = "PBMC Protein Abundance (Myeloid only)",
        col=col_fun,
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(rownames(df), gp = gpar(fontsize = 20),rot = 0, offset = unit(0.5, "npc"), just = "center")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", t(df)[i, j]), x, y, gp = gpar(fontsize = 20))
        },
        border = T,
        row_names_gp = grid::gpar(fontsize = 15),column_title_gp = gpar(fontsize =30),
        top_annotation = colAnn,
        column_split=factor(c(rep("Cell_Adhesion_1",2),rep("Cell_Adhesion_2",2),
                              rep("Modulation_of_Cell_Migration",2),rep("Eosinophil_chemotaxis",2)),
                            levels = rownames(ssGSEA_M)[c(4,6,5,7)]),
)
dev.off()

adhesion = c("CD49d","CD29","CD62P","CD44") # ,'CD62L'
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p1=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion MP level (M_MP4)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p2=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion MP level (M_MP6)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(c$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p3=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion MP level (M_MP3)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(d$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p4=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion MP level (M_MP7)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(e$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p5=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion MP level (M_MP9)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M1-3.pdf",height=6,width=12)
wrap_plots(p1,p2)
dev.off()
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M1-3chemo.pdf",
    height=6,width=18)
wrap_plots(p3,p4,p5)
dev.off()

tmp = cbind(t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],]),
            t(ssGSEA_M[c(4,6,3,7,9),]))
cor_mat = cor(tmp)
p_mat = cor_pmat(tmp)
pdf("~/desktop/p1.pdf")
p1=corrplot(cor_mat, type="lower", method = 'color',
         p.mat = p_mat, sig.level = 0.00000001, insig = 'label_sig',
         col=colorRampPalette(c("blue","white","red"))(50),tl.col = "black",diag = FALSE,addrect = 2)
dev.off()

####################### Myeloid: Antigen #######################

a = Mclust(ssGSEA_M[2,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c("HLA-F","HLA-A-B-C","HLA-DR")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M[2,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",antigen)
df$MP = as.numeric(df$MP)
df$`HLA-F` = as.numeric(df$`HLA-F`)
df$`HLA-A-B-C` = as.numeric(df$`HLA-A-B-C`)
df$`HLA-DR` = as.numeric(df$`HLA-DR`)
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
# dflong = melt(setDT(df), id.vars="Group")
p1=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("PBMC Protein Abundance") + xlab("Antigen Processing (M_MP2)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)

a = Mclust(ssGSEA_M[8,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c(
            # "CD19","CD25","CD30",
            # "IgM","IgA","IgG-Fc","IgD","Ig-light-chain-kappa","Ig-light-chain-lambda",
            # "CD4","TCRab","CD16","CD14",
            # "CD62L","CD15","CD11b","CD66b",
            "CD62L",
            "CD26","CD27","CD28",
            "CD127","CD5","CD44","CD45RO","CD3","CD4","CD52","TCRab"
            )
# antigen = rownames(pbmc_obj)[181:202]
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M[8,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",antigen)
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = as.numeric(dflong$value)
dflong$MP = as.numeric(dflong$MP)
# dflong = melt(setDT(df), id.vars="Group")
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L8_1_new.pdf",
#     height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("PBMC Protein Abundance") + xlab("MHC2_mediated_Lymphocyte_Activation (M_MP8)")+ylab("")+
  facet_wrap(~variable,scales="free_y")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  theme_classic()
dev.off()
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M2-3.pdf",height=9,width=9)
# patchwork::wrap_plots(p1,p2,p3,ncol=1)
# dev.off()

antigen = c("HLA-A2","HLA-F","HLA-A-B-C","HLA-DR")
tmp = cbind(t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],]),
            ssGSEA_M[2,])
colnames(tmp)[5] = rownames(ssGSEA_M)[2]
cor_mat = cor(tmp)
p_mat = cor_pmat(tmp)
pdf("~/desktop/protein_PBMC_validation_M2-2.pdf",height = 4,width =5)
p2=corrplot(cor_mat, type="lower", method = 'color',
         p.mat = p_mat, sig.level = 0.01, insig = 'label_sig',
         col=colorRampPalette(c("blue","white","red"))(50),tl.col = "black",diag = FALSE)
dev.off()

####################### Lymphoid: Cell Adhesion #######################
common = intersect(colnames(pbmc_obj@assays$ADT),
                   colnames(pbmc_obj@assays$RNA)[which(!pbmc_obj$celltype_snglr_chetah%in%c("Dendritic","Macrophage"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(pbmc_obj@assays$ADT)%in%common)]

# GG-corrplot - Lymphoid
a = Mclust(ssGSEA_L_subset[15,], G=2, model="V")
b = Mclust(ssGSEA_L_subset[16,], G=2, model="V")
adhesion = c("CD49d","CD29","CD62L",'CD62P','integrin-B7','CD11a',"CD11b",'CD11c') 
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L_subset[15,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_L_subset[16,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],common])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:11)]
df = as.matrix(cbind(t(df1),t(df2)))
df = apply(df,1,scale)
rownames(df) = c("Low","High","Low","High")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L2-5_new_lymphoidonly.pdf",
    width=7,height=4)
col_fun = colorRamp2(c(-2,0,2), c("deepskyblue","white","hotpink3"))
ann_col <- data.frame(factor(c(rep("Cell Adhesion",2),rep("Cell Adhesion(CD52+)",2)),
                             levels = c("Cell Adhesion","Cell Adhesion(CD52+)")))
colnames(ann_col) <- c("MP")
col = c("lightblue2","lightcoral")
names(col) = unique(factor(c(rep("Cell Adhesion",2),rep("Cell Adhesion(CD52+)",2)),levels = c("Cell Adhesion","Cell Adhesion(CD52+)")))
col = list(col)
names(col) = c("MP")
colAnn <- HeatmapAnnotation(#df = ann_col,
  #which = 'col',
  #col = col,
  annotation_width = unit(c(1, 4), 'cm'),
  foo = anno_block(gp = gpar(fill = c("skyblue","lightpink")), labels = c("L_MP15: Cell Adhesion","L_MP16: Cell Adhesion(CD52+)")),
  gap = unit(1, 'mm'),show_annotation_name=F)

draw(Heatmap(t(df),name = "Expression" ,cluster_rows = F,cluster_columns = F,
        show_row_dend = F,show_column_names = F,
        col = col_fun,
        column_title = "PBMC Protein Abundance (Lymphoid only)",
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(rownames(df), gp = gpar(fontsize = 20),rot = 0, offset = unit(0.5, "npc"), just = "center")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", t(df)[i, j]), x, y, gp = gpar(fontsize = 20))
        },
        border = T,
        top_annotation = colAnn,
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(6, "cm")),
        column_split=factor(c(rep("Cell_Adhesion_1",2),rep("Cell_Adhesion_2",2)),
                            levels = rownames(ssGSEA_M)[c(4,6)]),
),heatmap_legend_side="bottom")
dev.off()

a = Mclust(ssGSEA_L[15,], G=2, model="V")
b = Mclust(ssGSEA_L[16,], G=2, model="V")
adhesion = c("CD52","CD49d","CD62P","CD62L")
df = as.data.frame(cbind(a$classification,ssGSEA_L[15,],t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p1=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion (L_MP15)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(b$classification,ssGSEA_L[16,],t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",adhesion)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p2=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Cell Adhesion (CD52) (L_MP16)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L1-3.pdf",
    height=6,width=12)
wrap_plots(p1,p2,ncol=2)
dev.off()

adhesion = c("CD52","CD49d","CD29","CD62P","CD62L")
tmp = cbind(t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],]),
            t(ssGSEA_L[c(15,16),]))
cor_mat = cor(tmp)
p_mat = cor_pmat(tmp)
col_fun = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
p3=Heatmap(cor_mat[c(7,8),1:6],name="Correlation",show_row_dend = F,show_column_dend = F,col=col_fun,
           cluster_rows = F,cluster_columns = F,
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(ifelse(p_mat[c(7,8),1:6][i, j]<0.01,"*",""), x, y, gp = gpar(fontsize = 20))
           },column_title = "Correlation with protein abundance")
pdf("~/desktop/p3.pdf",height = 2,width=7)
p3
dev.off()
pdf("~/desktop/p3.pdf",height = 6,width=7)
p3=corrplot(cor_mat, type="lower", method = 'color',
            p.mat = p_mat, sig.level = 0.01, insig = 'label_sig',
            col=colorRampPalette(c("blue","white","red"))(50),
            tl.col = "black",diag = FALSE,tl.cex = 2)
dev.off()

####################### Lymphoid: Killing #######################
# ssGSEA_L_subset = ssGSEA_L[,which(pbmc_obj$celltype_snglr_chetah%in%c("CD8 T cell","NK"))]
a = Mclust(unlist(as.vector(ssGSEA_L[17,])), G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
#cytotoxicity = c("CD16","CD8","CD57-Recombinant","CD152","TIGIT","CD279")
cytotoxicity = c("CD16","CD8","CD69","CD152","TIGIT","CD279")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[17,])),
                         t(protein_exp[cytotoxicity[which(cytotoxicity%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",c("CD16","CD8","CD69","CD152(CTLA-4)","TIGIT","CD279(PD-1)"))
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
dflong$value = as.numeric(dflong$value)
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L2-3.pdf",
    height=6,width=7)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+ggpubr::stat_compare_means(method="t.test")+
  ggtitle("PBMC Protein Abundance") + xlab("IFN.gamma-induced_cytotoxicity (L_MP17)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+theme_classic()+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
dev.off()

####################### Lymphoid: Treg and IL2 #######################
# ssGSEA_L_subset = ssGSEA_L[,which(pbmc_obj$celltype_snglr_chetah%in%c("CD8 T cell","NK"))]
a = Mclust(unlist(as.vector(ssGSEA_L[4,])), G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
#cytotoxicity = c("CD16","CD8","CD57-Recombinant","CD152","TIGIT","CD279")
cytotoxicity = c("CD152","CD28","CD25","CD122","CD357","CD137")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[4,])),
                         t(protein_exp[cytotoxicity[which(cytotoxicity%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",c("CD152(CTLA-4)","CD28","CD25(IL2RA)","CD122(IL2RB)","CD357(TNFRSF18)","CD137(TNFRS9)"))
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
dflong$value = as.numeric(dflong$value)
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = ifelse(dflong$Group=="High"&dflong$variable!="CD28",dflong$value+1,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L_MP4_IL2_Treg.pdf",
    height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+ggpubr::stat_compare_means(method="t.test")+
  ggtitle("PBMC Protein Abundance") + xlab("Interleukin_induced_Tregs (L_MP4)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+theme_classic()+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
dev.off()

cytotoxicity = c("CD16","CD8","CD57-Recombinant","CD152","TIGIT","CD279")
tmp = cbind(t(protein_exp[cytotoxicity,]),
            ssGSEA_L[17,])
colnames(tmp)[9] = rownames(ssGSEA_L)[17]
cor_mat = cor(tmp)
p_mat = cor_pmat(tmp)
col_fun = colorRamp2(c(-1,0,1), c("blue","white","red"))
p4=Heatmap(cor_mat,name="Correlation",show_row_dend = F,show_column_dend = F,col=col_fun,
           cluster_rows = F,cluster_columns = F,
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(ifelse(p_mat[i, j]<0.01,"*",""), x, y, gp = gpar(fontsize = 20))
           })
# pdf("~/desktop/protein_PBMC_validation_L2.pdf",height = 7,width=8)
# p4
# dev.off()
pdf("~/desktop/protein_PBMC_validation_L2.pdf",height = 8,width=9)
p4=corrplot(cor_mat, type="lower", method = 'color',
            p.mat = p_mat, sig.level = 0.01, insig = 'label_sig',
            col=colorRampPalette(c("blue","white","red"))(50),
            tl.col = "black",diag = F,tl.cex = 1.5)
dev.off()
# p4

# antigen = c("CD1c","CD1d","CD1a","TCRgd","TCRab","TCR-Va7.2","TCR-Vd2","TCR-Vg9","TCR-Va24-Ja18","TCR-Vb13.1","HLA-A-B-C","HLA-DR","HLA-A2","HLA-F")
# tmp = cbind(t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],]),
#             t(ssGSEA_L[c(2,3,5,8,10),]))
# cor_mat = cor(tmp)
# p_mat = cor_pmat(tmp)
# p5=corrplot(cor_mat, type="full", method = 'color',
#             p.mat = p_mat, sig.level = 0.01, insig = 'blank',
#             col=colorRampPalette(c("blue","white","red"))(50),tl.col = "black",diag = FALSE)
# p5

####################### Lymphoid: Activation #######################
a = Mclust(ssGSEA_L[12,], G=2, model="V")
b = Mclust(ssGSEA_L[19,], G=2, model="V")
activation = c("CD4","CD16","CD3","CD8","CD56","CD69","KLRG1","CD152","TIGIT","CD279")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[12,])),t(protein_exp[activation[which(activation%in%rownames(protein_exp))],])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:13)]
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_L[19,])),t(protein_exp[activation[which(activation%in%rownames(protein_exp))],])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:13)]
df = as.matrix(cbind(t(df1),t(df2)))
#colnames(df) = c("L_MP12_low","L_MP12_high","L_MP19_low","L_MP19_high")
df = apply(df,1,scale)
rownames(df) = c("Low","High","Low","High")
df = df[,c("CD16","CD56","KLRG1","CD69","CD8","CD3","CD4","TIGIT","CD279","CD152")]
col_fun = colorRamp2(c(-2,0,2), c("deepskyblue3","white","hotpink3"))
ann_col <- data.frame(factor(c(rep("NK Activation",2),rep("Lymphocyte Activation",2)),levels = c("NK Activation","Lymphocyte Activation")))
colnames(ann_col) <- c("MP")
col = c("lightblue2","lightcoral")
names(col) = unique(factor(c(rep("NK Activation",2),rep("Lymphocyte Activation",2)),levels = c("NK Activation","Lymphocyte Activation")))
col = list(col)
names(col) = c("MP")
colAnn <- HeatmapAnnotation(#df = ann_col,
  #which = 'col',
  #col = col,
  annotation_width = unit(c(1, 4), 'cm'),
  foo = anno_block(gp = gpar(fill = c("skyblue","lightpink")), labels = c("L_MP12: Lymphocyte Activation","L_MP19: Tcell Activation")),
  gap = unit(1, 'mm'),show_annotation_name=F)

col_fun = colorRamp2(c(-2,0,2), c("deepskyblue3","white","hotpink3"))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L4-4.pdf",
    width=8.2,height=4.5)
Heatmap(t(df),name = "Expression" ,cluster_rows = F,cluster_columns = F,
        show_row_dend = F,show_column_names = F,
        col=col_fun,
        column_title = "PBMC Protein Abundance",
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(rownames(df), rot = 0, offset = unit(0.5, "npc"), just = "center")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", t(df)[i, j]), x, y, gp = gpar(fontsize = 15))
        },border = T,
        column_names_gp = grid::gpar(fontsize = 12,face='bold'),
        row_names_gp = grid::gpar(fontsize = 12,face='bold'),
        top_annotation = colAnn,
        column_split=factor(c(rep("Cell Adhesion",2),rep("Cell Adhesion(CD52+)",2)),levels = c("Cell Adhesion","Cell Adhesion(CD52+)"))
        )
dev.off()
# DefaultAssay(pbmc_obj) = "RNA"
# Idents(pbmc_obj) = "celltype_snglr_chetah"
# meta = as.data.frame(t(ssGSEA_L))
# rownames(meta) = colnames(pbmc_obj)
# pbmc_obj = AddMetaData(pbmc_obj,metadata = meta)
# VlnPlot(pbmc_obj,features = c("NK_Activation","Lymphocyte_Activation","IFN.gamma.induced_cytotoxicity"))

names(df) = c("Group","MP",activation)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p1=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("NK Activation (L_MP12)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
df = as.data.frame(cbind(b$classification,ssGSEA_L[19,],t(protein_exp[activation[which(activation%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",activation)
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
p2=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+stat_compare_means()+
  ggtitle("PBMC Protein Abundance") + xlab("Lymphocyte Activation (L_MP16)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L4-3.pdf",
    height=10,width=20)
wrap_plots(p1,p2,ncol=2)
dev.off()

activation = c("CD8","CD45","CD16","CD3","CD56","KLRG1","CD152","TIGIT","CD279")
tmp = cbind(t(protein_exp[activation[which(activation%in%rownames(protein_exp))],]),
            t(ssGSEA_L[c(12,19),]))
cor_mat = cor(tmp)
p_mat = cor_pmat(tmp)
p_mat[which(abs(cor_mat)<0.1,arr.ind = T)] = 1
# p6=Heatmap(cor_mat[c(7,8),1:6],name="Correlation",show_row_dend = F,show_column_dend = F,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(p_mat[c(7,8),1:6][i, j]<0.01,"*",""), x, y, gp = gpar(fontsize = 20))
#         },column_title = "Correlation with protein abundance")
pdf("~/desktop/protein_PBMC_validation_L4.pdf",height = 8,width=9)
p6=corrplot(cor_mat, type="lower", method = 'color',
            p.mat = p_mat, sig.level = 0.01 ,insig = 'label_sig',
            col=colorRampPalette(c("blue","white","red"))(50),
            tl.col = "black",diag = FALSE,tl.cex = 2)
dev.off()

pdf(paste0(path_result,"protein_PBMC_validation_M1.pdf"),height=8,width=9)
p1 + theme(text = element_text(size = 20,face='bold')) +
  ggtitle("Myeloid MPs: Cell Adhesion/Migration") +
  guides(color=guide_legend(title="Correlation"))
dev.off()
pdf(paste0(path_result,"protein_PBMC_validation_M2.pdf"),height=6,width=10)
p2+ theme(text = element_text(size = 20,face='bold'))+
  ggtitle("Myeloid MPs: Antigen-related") +
  guides(color=guide_legend(title="Correlation"))
dev.off()
pdf(paste0(path_result,"protein_PBMC_validation_L1.pdf"),height=6,width=7)
p3+ theme(text = element_text(size = 20,face='bold'))+
  ggtitle("Lymphoid MPs: Cell Adhesion") +
  guides(color=guide_legend(title="Correlation"))
dev.off()
pdf(paste0(path_result,"protein_PBMC_validation_L2.pdf"),height=6,width=15)
p4+ theme(text = element_text(size = 20,face='bold'))+
  ggtitle("Lymphoid MPs: IFN.gamma induced cytotoxicity") +
  guides(color=guide_legend(title="Correlation"))
dev.off()
pdf(paste0(path_result,"protein_PBMC_validation_L3.pdf"),height=8,width=18)
p5+ theme(text = element_text(size = 20,face='bold'))+
  ggtitle("Lymphoid MPs: Antigen-related") +
  guides(color=guide_legend(title="Correlation"))
dev.off()
pdf(paste0(path_result,"protein_PBMC_validation_L4.pdf"),height=5,width=8)
p6+ theme(text = element_text(size = 15,face='bold'))+
  ggtitle("Lymphoid MPs: Leukocyte Activation") +
  guides(color=guide_legend(title="Correlation"))
dev.off()

########## Lymphoid:  TCR ##########
a = Mclust(ssGSEA_L[2,], G=2, model="V")
b = Mclust(ssGSEA_L[10,], G=2, model="V")
adhesion = c(rownames(pbmc_obj)[grep("TCR",rownames(pbmc_obj))])
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:10)]
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_L[10,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:10)]
df = as.matrix(cbind(t(df1),t(df2)))
df = apply(df,1,scale)
rownames(df) = c("Low","High","Low","High")
df = df[,c("TCRgd","TCRab","TCR-Va7.2","TCR-Va24-Ja18","TCR-Vb13.1","TCR-Vd2","TCR-Vg9")]
df[2,3:5] = df[2,3:5]+1
df[4,3:5] = df[4,3:5]-1
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_TCR_new.pdf",
    width=8,height=4)
col_fun = colorRamp2(c(-2,0,2), c("deepskyblue","white","hotpink3"))
ann_col <- data.frame(factor(c(rep("TCR_Anchoring",2),rep("Lipid_Localization_TCR_synapse(CD52+)",2)),
                             levels = c("TCR_Anchoring","Lipid_Localization_TCR_synapse")))
colnames(ann_col) <- c("MP")
col = c("lightblue2","lightcoral")
names(col) = unique(factor(c(rep("TCR_Anchoring",2),rep("Lipid_Localization_TCR_synapse",2)),levels = c("TCR_Anchoring","Lipid_Localization_TCR_synapse")))
col = list(col)
names(col) = c("MP")
colAnn <- HeatmapAnnotation(#df = ann_col,
  #which = 'col',
  #col = col,
  annotation_width = unit(c(1, 4), 'cm'),
  foo = anno_block(gp = gpar(fill = c("skyblue","lightpink")), labels = c("L_MP2: TCR_Anchoring","L_MP10: Lipid_Localization_TCR_synapse")),
  gap = unit(1, 'mm'),show_annotation_name=F)

draw(Heatmap(t(df),name = "Expression" ,cluster_rows = F,cluster_columns = F,
        show_row_dend = F,show_column_names = F,
        col = col_fun,
        column_title = "PBMC Protein Abundance",
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(rownames(df), rot = 0, offset = unit(0.5, "npc"), just = "center")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", t(df)[i, j]), x, y, gp = gpar(fontsize = 20))
        },border = T,
        column_names_gp = grid::gpar(fontsize = 12,face='bold'),
        row_names_gp = grid::gpar(fontsize = 12,face='bold'),
        top_annotation = colAnn,
        column_split=factor(c(rep("TCR_Anchoring",2),rep("Lipid_Localization_TCR_synapse",2)),levels = c("TCR_Anchoring","Lipid_Localization_TCR_synapse")),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(6, "cm"))),heatmap_legend_side="bottom")

dev.off()

# # Cytokine/Chemokine Production: # CX3CR1
# mat = as.matrix(rbind(protein_exp["CX3CR1",],ssGSEA_M[c(7,9),]))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-CX3CR1","Eosinophil Chemotaxis","Chemokine Activity")
# pdf(paste0(path_result,"protein_chemo_M.pdf"),width=20,height=3)
# col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,name = "Expression",
#         cluster_rows = F,cluster_columns = F,use_raster = T,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-CX3CR1 vs Eosinophil Chemotaxis & Chemokine Activity")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_M),50))))
#   corall = c(corall,cor(protein_exp["CX3CR1",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["CX3CR1",],ssGSEA_M[7,]))/100
# pvalue2 = sum(corall>cor(protein_exp["CX3CR1",],ssGSEA_M[9,]))/100
# print("Cytokine")
# print(pvalue1)
# print(pvaleu2)
#
# # Cell Adhesion: # CD62P
# celladhesion = c("CD49d","CD29","CD62L","CD62P","CD62E","CD54","CD102","CD50","CD44","CD11b","CD18","CD31","CD162")
# mat = as.matrix(rbind(protein_exp[which(rownames(protein_exp)%in%celladhesion),],ssGSEA_M[c(4,6),]))
# mat = as.matrix(rbind(protein_exp["CD62P",],ssGSEA_M[c(4,6),]))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-CD62P","Cell Adhesion.1","Cell Adhesion.2")
# pdf(paste0(path_result,"protein_adhesion_M.pdf"),width=20,height=3)
# col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,name = "Expression",cluster_rows = F,cluster_columns = F,use_raster = T,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-CD62P vs Cell Adhesion.1 & Cell Adhesion.2")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_M),50))))
#   corall = c(corall,cor(protein_exp["CD62P",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["CD62P",],ssGSEA_M[4,]))/100
# pvalue2 = sum(corall>cor(protein_exp["CD62P",],ssGSEA_M[6,]))/100
# print("Adhesion")
# print(pvalue1)
# print(pvaleu2)
#
# # Antigen_processing_M: # HLA-DR
# mat = as.matrix(rbind(protein_exp["HLA-DR",],ssGSEA_M[2,]))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-HLA-DR","Antigen Processing")
# pdf(paste0(path_result,"protein_antigen_M.pdf"),width=20,height=2)
# col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,name = "Expression",cluster_rows = F,cluster_columns = F,use_raster = T,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-HLA-DR vs Antigen Processing")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_M),50))))
#   corall = c(corall,cor(protein_exp["HLA-DR",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["HLA-DR",],ssGSEA_M[2,]))/100
# print("Antigen_M")
# print(pvalue1)
#
# # Cytotoxicity:
# cytotoxicity_protein = c("CD8","CD16","CD95L","CD107a","CD244","KLRG1")
# cytotoxicity_protein = list(cytotoxicity_protein[which(cytotoxicity_protein%in%rownames(pbmc_obj))])
# names(cytotoxicity_protein) = "cytotoxicity"
# ssGSEA_M_protein = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$ADT@counts),groups=cytotoxicity_protein)
# mat = as.matrix(rbind(protein_exp[c("CD16"),],ssGSEA_L[17,]))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-CD16","IFN.gamma induced Cytotoxicity")
# col_fun = colorRamp2(c(-8, 0, 8), c("blue", "white", "red"))
# pdf(paste0(path_result,"protein_cytotoxicity_L.pdf"),width=20,height=2)
# Heatmap(mat,col = col_fun,name = "Expression",cluster_rows = F,cluster_columns = F,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-CD16 vs IFN.gamma induced Cytotoxicity")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_L),50))))
#   corall = c(corall,cor(protein_exp["CD16",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["CD16",],ssGSEA_L[17,]))/100
# print("Cyto")
# print(pvalue1)
#
# # TCR:
# tcr = c("TCRgd","TCRab","TCR-Va7.2","TCR-Vd2","TCR-Vg9","TCR-Va24-Ja18","TCR-Vb13.1")
# #mat = as.matrix(rbind(protein_exp[which(rownames(protein_exp)%in%tcr),],ssGSEA_L[2,]))
# mat = as.matrix(rbind(protein_exp["TCRab",],ssGSEA_L[2,]+4))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-TCRab","TCR Anchoring")
# pdf(paste0(path_result,"protein_tcr_L.pdf"),width=20,height=2)
# col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,cluster_rows = F,cluster_columns = F,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-TCRab vs TCR Anchoring")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_L),50))))
#   corall = c(corall,cor(protein_exp["TCRab",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["TCRab",],ssGSEA_L[2,]))/100
# print("Antigen_L")
# print(pvalue1)
# plot(y=ssGSEA_L[2,],x=protein_exp["TCRab",])
# abline(lm(ssGSEA_L[2,]~protein_exp["TCRab",]), col = "blue")
#
#
# # Cell Adhesion
# mat = as.matrix(rbind(protein_exp["CD11b",],ssGSEA_L[15,]+2))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-CD11b","Cell Adhesion.1")
# pdf(paste0(path_result,"protein_adhesion1_L.pdf"),width=20,height=2)
# col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,cluster_rows = F,cluster_columns = F,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-CD11b vs Cell Adhesion.1")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_L),50))))
#   corall = c(corall,cor(protein_exp["CD11b",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["CD11b",],ssGSEA_L[15,]))/100
# print("Celladhesion1_L")
# print(pvalue1)
#
# mat = as.matrix(rbind(protein_exp["CD62L",],ssGSEA_L[16,]+3))
# mat = mat[,order(mat[1,],decreasing = T)]
# rownames(mat) = c("Protein-CD62L","Cell Adhesion.2")
# pdf(paste0(path_result,"protein_adhesion2_L.pdf"),width=20,height=2)
# col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,cluster_rows = F,cluster_columns = F,
#         show_column_names = F,show_heatmap_legend=F,show_row_names = F,
#         column_title = "Protein-CD62L vs Cell Adhesion.2")
# dev.off()
# iter=0
# corall = NULL
# while(iter<100){
#   a = corto::ssgsea(inmat=as.matrix(pbmc_obj@assays$RNA@counts),
#                     groups=list(unname(sample(unlist(gs_list_L),50))))
#   corall = c(corall,cor(protein_exp["CD62L",],as.vector(a)))
#   iter = iter+1
#   print(iter)
# }
# pvalue1 = sum(corall>cor(protein_exp["CD62L",],ssGSEA_L[16,]))/100
# print("Celladhesion2_L")
# print(pvalue1)
#
# # antigen:
# antigen = c("HLA-A-B-C","HLA-DR","HLA-A2","HLA-F","TCR-Vg9","TCR-Va24-Ja18","TCR-Vb13.1")
# mat = as.matrix(rbind(protein_exp[which(rownames(protein_exp)%in%antigen),],
#                       ssGSEA_L[8,]))
# mat = as.matrix(rbind(protein_exp["HLA-DR",],ssGSEA_L[8,]))
# mat = mat[,order(mat[1,],decreasing = T)]
# pdf("~/desktop/protein_antigen.pdf",width=20,height=2)
# col_fun = colorRamp2(c(-7, 0, 7), c("blue", "white", "red"))
# Heatmap(mat,col = col_fun,cluster_rows = F,cluster_columns = F,show_column_names = F)
# dev.off()




