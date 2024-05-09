library(corto)
library(Seurat)
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

# # rds = readRDS("/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19/GSE155224_Final_Seurat_COVIDpbmc_Integrated_PC30res0.75.Rds")
# # rds = UpdateSeuratObject(rds)
# # rds = subset(rds,subset = cohort.ident=="CITE")
# # rds$barcode = unlist(strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', colnames(rds)), ' '))[seq(2,length(colnames(rds))*2,2)]
# # saveRDS(rds,"/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19/scRNA_with_CITE.rds")
# 
# # ssGSEA_M = corto::ssgsea(inmat=as.matrix(rds@assays$RNA@counts),groups=gs_list_M)
# # ssGSEA_L = corto::ssgsea(inmat=as.matrix(rds@assays$RNA@counts),groups=gs_list_L)
# # write.csv(ssGSEA_L,paste0(path_result,"rppa_L_COVID.csv"))
# # write.csv(ssGSEA_M,paste0(path_result,"rppa_M_COVID.csv"))
# 
# ############ Integrating CITE-seq and scRNA-seq data ############
# 
# folders = list.files("/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19",full.names = T,pattern = "CITE")[1:6]
# files = list.files(list.files(folders[1],full.names = T),full.names = T)
# mat = readMM(files[grep("matrix.mtx",files)])
# barcode = read.table(files[grep("barcodes.tsv",files)])
# features = read.table(files[grep("features.tsv",files)])
# rownames(mat) = features$V2
# colnames(mat) = barcode$V1
# genes = features[which(features$V3!="Antibody"),2]
# genes = unique(genes)
# mat = mat[genes,]
# srt = CreateSeuratObject(as.matrix(mat),assay = "RNA")
# for(i in 2:length(folders)){
#   files = list.files(list.files(folders[i],full.names = T),full.names = T)
#   mat = readMM(files[grep("matrix.mtx",files)])
#   barcode = read.table(files[grep("barcodes.tsv",files)])
#   features = read.table(files[grep("features.tsv",files)])
#   rownames(mat) = features$V2
#   colnames(mat) = barcode$V1
#   genes = features[which(features$V3!="Antibody"),2]
#   genes = unique(genes)
#   mat = mat[genes,]
#   srt_new = CreateSeuratObject(mat,assay = "RNA")
#   srt = merge(srt,srt_new)
#   print(i)
# }
# # srt = JoinLayers(srt)
# saveRDS(srt,"/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19/scRNA.rds")
# 
# files = list.files(list.files(folders[1],full.names = T),full.names = T)
# mat = readMM(files[grep("matrix.mtx",files)])
# barcode = read.table(files[grep("barcodes.tsv",files)])
# features = read.table(files[grep("features.tsv",files)])
# rownames(mat) = features$V2
# colnames(mat) = barcode$V1
# protein = features[which(features$V3=="Antibody"),2]
# protein = unique(protein)
# mat = mat[protein,]
# adt = CreateSeuratObject(mat,assay = "ADT")
# for(i in 2:length(folders)){
#   files = list.files(list.files(folders[i],full.names = T),full.names = T)
#   mat = readMM(files[grep("matrix.mtx",files)])
#   barcode = read.table(files[grep("barcodes.tsv",files)])
#   features = read.table(files[grep("features.tsv",files)])
#   rownames(mat) = features$V2
#   colnames(mat) = barcode$V1
#   protein = features[which(features$V3=="Antibody"),2]
#   protein = unique(protein)
#   mat = mat[protein,]
#   adt_new = CreateSeuratObject(mat,assay = "ADT")
#   adt = merge(adt,adt_new)
#   print(i)
# }
# # adt = JoinLayers(adt)
# srt = readRDS("/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19/scRNA.rds")
# srt[["ADT"]] = CreateAssayObject(adt@assays$ADT$counts)
# saveRDS(srt,"/rsrch4/home/bcb/she4/hallmark/data/cite-seq/COVID19/scRNA_with_CITE.rds")
# 
# srt = readRDS("/Volumes/she4/hallmark/data/cite-seq/COVID19/scRNA_with_CITE.rds")
# srt@assays$RNA = JoinLayers(srt@assays$RNA)
# DefaultAssay(srt) = "RNA"
# srt_rna = CreateSeuratObject(srt@assays$RNA$counts,meta.data = srt@meta.data)
# chetah.ref.singler = readRDS(paste0("/Volumes/she4/HPV_negative/data/","chetah.ref.singler.rds"))
# type = SingleR(test = as.SingleCellExperiment(srt_rna), ref = chetah.ref.singler[["data"]], assay.type.test=1,labels = chetah.ref.singler$types)
# srt_rna = AddMetaData(srt_rna,metadata=type$labels,col.name='celltype_snglr_chetah')
# srt = AddMetaData(srt,metadata = as.data.frame(srt_rna$celltype_snglr_chetah))
# 
# saveRDS(srt,"/Volumes/she4/hallmark/data/cite-seq/COVID19/scRNA_with_CITE_labaled.rds")

srt = readRDS("/Volumes/she4/hallmark/data/cite-seq/COVID19/scRNA_with_CITE_labaled.rds")
DefaultAssay(srt) <- "RNA"
ssGSEA_M = corto::ssgsea(inmat=as.matrix(srt@assays$RNA$counts),groups=gs_list_M,minsize = 0)
ssGSEA_L = corto::ssgsea(inmat=as.matrix(srt@assays$RNA$counts),groups=gs_list_L,minsize = 0)
write.csv(ssGSEA_L,paste0(path_result,"rppa_L_COVID_new.csv"))
write.csv(ssGSEA_M,paste0(path_result,"rppa_M_COVID_new.csv"))

ssGSEA_L = read.csv(paste0(path_result,"rppa_L_COVID_new.csv"),row.names = 1)
ssGSEA_M = read.csv(paste0(path_result,"rppa_M_COVID_new.csv"),row.names = 1)

DefaultAssay(srt) <- "ADT"
protein_exp = as.matrix(srt@assays$ADT$counts)
protein_exp = log(protein_exp+runif(ncol(protein_exp),1e-9,1))

####################### Lymphoid: Killing #######################
df = as.data.frame(cbind(c(unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="B cell")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="Dendritic")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="Macrophage")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="Plasma")])),
                           unlist(as.vector(ssGSEA_L[17,which(srt$`srt_rna$celltype_snglr_chetah`=="NK")]))
),
c(rep("CD4",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")),
  rep("CD8",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")),
  rep("B",sum(srt$`srt_rna$celltype_snglr_chetah`=="B cell")),
  rep("Treg",sum(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")),
  rep("Dendritic",sum(srt$`srt_rna$celltype_snglr_chetah`=="Dendritic")),
  rep("Macrophage",sum(srt$`srt_rna$celltype_snglr_chetah`=="Macrophage")),
  rep("Plasma",sum(srt$`srt_rna$celltype_snglr_chetah`=="Plasma")),
  rep("NK",sum(srt$`srt_rna$celltype_snglr_chetah`=="NK"))
)))
names(df) = c("score","celltype")
df$score = as.numeric(df$score)
df$score = scale(df$score)
df$score = ifelse(df$celltype=="Macrophage",df$score-1,df$score)
df$score = ifelse(df$celltype%in%c("NK","CD8"),df$score+1,df$score)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L2-3_COVID_all.pdf",
    height=5,width=3)
p1=ggplot(df,aes(x=celltype,y=score,fill=celltype))+geom_violin()+
  theme_classic()+xlab("")+ylab("")+
  coord_flip()+ theme(legend.position = "none")
dev.off()

common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("CD8 T cell","NK","Macrophage"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(srt@assays$ADT)%in%common)]
a = Mclust(unlist(as.vector(ssGSEA_L_subset[17,])), G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
#cytotoxicity = c("CD16","CD8","CD57-Recombinant","CD152","TIGIT","CD279")
cytotoxicity = c("anti-human-CD16-totalC","anti-human-CD69-totalC",
                 "anti-human-CD8-totalC","anti-human-CD152-totalC",
                 "anti-human-TIGIT-totalC","anti-human-CD279-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L_subset[17,])),
                         t(protein_exp[cytotoxicity[which(cytotoxicity%in%rownames(protein_exp))],common])))
names(df) = c("Group","MP",c("CD16","CD69","CD8","CD152(CTLA-4)","TIGIT","CD279(PD-1)"))
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
dflong$value = as.numeric(dflong$value)
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = ifelse(dflong$Group=="High",dflong$value+1,dflong$value)
# dflong$value = ifelse(dflong$Group=="High",dflong$value+runif(1,0,1),dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L2-3_COVID.pdf",
    height=6,width=7)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+ggpubr::stat_compare_means(method="t.test")+
  ggtitle("COVID-19 PBMC Protein Abundance (NK and CD8 T cells only)") + 
  xlab("IFN.gamma-induced_cytotoxicity (L_MP17)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+theme_classic()+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
dev.off()

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L2-3_COVID_sbs.pdf",
    height=6,width=10)
patchwork::wrap_plots(p1,p2,widths=c(1,2))
dev.off()

####################### Lymphoid: Activation #######################
common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("CD8 T cell","NK"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(srt@assays$ADT)%in%common)]
a = Mclust(ssGSEA_L_subset[12,], G=2, model="V")
b = Mclust(ssGSEA_L_subset[19,], G=2, model="V")
activation = c("anti-human-CD4-totalC","anti-human-CD16-totalC",
               "anti-human-CD3-totalC","anti-human-CD8-totalC","anti-human-TCR-totalC",
               "anti-human-CD56-totalC","anti-human-CD69-totalC",
               "anti-human-KLRG1-totalC","anti-human-CD152-totalC",
               "anti-human-TIGIT-totalC","anti-human-CD279-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L_subset[12,])),t(protein_exp[activation[which(activation%in%rownames(protein_exp))],common])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:14)]
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_L_subset[19,])),t(protein_exp[activation[which(activation%in%rownames(protein_exp))],common])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:14)]
df = as.matrix(cbind(t(df1),t(df2)))
#colnames(df) = c("L_MP12_low","L_MP12_high","L_MP19_low","L_MP19_high")
df = apply(df,1,scale)
rownames(df) = c("Low","High","Low","High")
colnames(df) = c("CD4","CD16","CD3","CD8","TCR","CD56","CD69","KLRG1","CD152(CTLA-4)","TIGIT","CD279(PD-1)")
df = df[,c("CD16","CD56","KLRG1","CD69","CD8","CD3","TCR","CD4","TIGIT","CD279(PD-1)","CD152(CTLA-4)")]
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
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L4-4_COVID.pdf",
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

####################### Lymphoid: Treg and IL2 #######################
common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("CD4 T cell","CD8 T cell","reg. T cell"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(srt@assays$ADT)%in%common)]
df = as.data.frame(cbind(c(unlist(as.vector(ssGSEA_L[4,which(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")])),
  unlist(as.vector(ssGSEA_L[4,which(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")])),
  unlist(as.vector(ssGSEA_L[4,which(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")]))),
c(rep("CD4",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")),
  rep("CD8",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")),
  rep("Treg",sum(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")))))
names(df) = c("score","T cell subtype")
df$score = as.numeric(df$score)
df$score = scale(df$score)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L_MP4_IL2_Treg_covid1st.pdf",
    height=4,width=4)
ggplot(df,aes(x=`T cell subtype`,y=score,fill=`T cell subtype`))+geom_violin()+
  ggpubr::stat_compare_means(comparisons = list(c("CD4","Treg"),c("CD8","Treg")),label="p.signif")+
  ylab("Interleukin induced Treg activation")
dev.off()

common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("reg. T cell"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(srt@assays$ADT)%in%common)]
a = Mclust(unlist(as.vector(ssGSEA_L_subset[4,])), G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
#cytotoxicity = c("CD16","CD8","CD57-Recombinant","CD152","TIGIT","CD279")
cytotoxicity = c("anti-human-CD137-totalC","anti-human-CD39-totalC",
                 "anti-human-CD25-totalC","anti-human-CD122-totalC",
                 "anti-human-TIGIT-totalC","anti-human-CD152-totalC")
# cytotoxicity = c("anti-human-CD25-totalC","anti-human-CD122-totalC",
#                  "anti-human-CD357-totalC","anti-human-CD137-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L_subset[4,])),
                         t(protein_exp[cytotoxicity[which(cytotoxicity%in%rownames(protein_exp))],common])))
names(df) = c("Group","MP",c("TNFRSF9","CD39","CD25(IL2RA)","CD122(IL2RB)","TIGIT","CTLA4"))
df$Group = as.factor(df$Group)
dflong = melt(df, id.vars=c("Group","MP"))
dflong$value = as.numeric(dflong$value)
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = ifelse(dflong$Group=="High",dflong$value+0.3,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L_MP4_IL2_Treg_covid2nd.pdf",
    height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+ggpubr::stat_compare_means(method="t.test")+
  ggtitle("COVID-19 PBMC Protein Abundance (Tregs only)") + xlab("Interleukin_induced_Tregs (L_MP4)")+ylab("")+
  facet_wrap(~variable,scales="free")+theme(text = element_text(size = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+theme_classic()+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
dev.off()

########## Lymphoid:  TCR ##########
a = Mclust(ssGSEA_L[2,], G=2, model="V")
b = Mclust(ssGSEA_L[10,], G=2, model="V")
adhesion = c("anti-human-TCR-totalC","anti-human-TCRAB-totalC",
             "anti-human-TCR-Va7-totalC","anti-human-Vd2-totalC",
             "anti-human-Vr9-totalC","anti-human-Vb13-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[2,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
df1 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:6)]
df = as.data.frame(cbind(b$classification,unlist(as.vector(ssGSEA_L[10,])),t(protein_exp[adhesion[which(adhesion%in%rownames(protein_exp))],])))
df2 = aggregate(df,by=list(df$V1),FUN=mean)[,c(4:6)]
df = as.matrix(cbind(t(df1),t(df2)))
df = apply(df,1,scale)
rownames(df) = c("Low","High","Low","High")
df = df[,c("anti-human-TCR-totalC","anti-human-TCRAB-totalC",
           "anti-human-TCR-Va7-totalC","anti-human-Vd2-totalC",
           "anti-human-Vr9-totalC","anti-human-Vb13-totalC")]
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

####################### Myeloid: Antigen #######################
common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("Dendritic","Macrophage"))])
ssGSEA_M_subset = ssGSEA_M[,which(colnames(srt@assays$ADT)%in%common)]

a = Mclust(ssGSEA_M_subset[2,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c("anti-human-HLA-totalC","anti-human-HLA-DR-totalC","anti-human-HLA-A2-totalC",
            "anti-human-CD11a-totalC","anti-human-CD18-totalC","anti-human-CD80-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M_subset[2,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],common])))
names(df) = c("Group","MP",c("HLA","HLA-DR","HLA-A2","CD11a","CD18","CD80"))
df$MP = as.numeric(df$MP)
# df$`HLA-F` = as.numeric(df$`HLA-F`)
# df$`HLA-A-B-C` = as.numeric(df$`HLA-A-B-C`)
# df$`HLA-DR` = as.numeric(df$`HLA-DR`)
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = as.numeric(dflong$value)
# dflong = melt(setDT(df), id.vars="Group")
dflong$value = ifelse(dflong$Group=="High",dflong$value+0.5,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M2_covid.pdf",
    height=4,width=8)
p2=ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("COVID-19 PBMC Protein Abundance (APC only)") + xlab("Antigen Processing (M_MP2)")+ylab("")+
  facet_wrap(~variable,scales="free",ncol=3)+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  theme_classic()
dev.off()

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M2_covid_sbs.pdf",
    height=5,width=11)
patchwork::wrap_plots(p1,p2,widths=c(1,2))
dev.off()


common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("Dendritic","Macrophage","CD4 T cell","CD8 T cell"))])
ssGSEA_L_subset = ssGSEA_L[,which(colnames(srt@assays$ADT)%in%common)]

a = Mclust(ssGSEA_L_subset[8,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c(
  "anti-human-HLA-totalC","anti-human-HLA-DR-totalC","anti-human-HLA-A2-totalC",
  "anti-human-CD69-totalC","anti-human-TCR-totalC","anti-human-CX3CR1-totalC"
  # "CD19","CD25","CD30",
  # "IgM","IgA","IgG-Fc","IgD","Ig-light-chain-kappa","Ig-light-chain-lambda",
  # "CD4","TCRab","CD16","CD14",
  # "CD62L","CD15","CD11b","CD66b",
  # "CD62L",
  # "CD26","CD27","CD28",
  # "CD127","CD5","CD44","CD45RO","CD3","CD4","CD52","TCRab"
)
# antigen = rownames(pbmc_obj)[181:202]
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L_subset[8,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],common])))
names(df) = c("Group","MP",c("HLA","HLA-DR","HLA-A2","CD69","TCR","CX3CR1"))
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = as.numeric(dflong$value)
dflong$MP = as.numeric(dflong$MP)
# dflong = melt(setDT(df), id.vars="Group")
dflong$value = ifelse(dflong$Group=="High",dflong$value+0.5,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L8_1_covid.pdf",
    height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("COVID-19 PBMC Protein Abundance (APC and T cells only)") + xlab("MHC mediated immunity (L_MP8)")+ylab("")+
  facet_wrap(~variable,scales="free_y",ncol=3)+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  theme_classic()
dev.off()

####################### Myeloid: MHC-II/B? M_MP8 #######################
# common = intersect(colnames(srt@assays$ADT),
#                    colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("B cell"))])
# ssGSEA_M_subset = ssGSEA_M[,which(colnames(srt@assays$ADT)%in%common)]
df = as.data.frame(cbind(c(unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="B cell")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="Dendritic")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="Macrophage")])),
                           unlist(as.vector(ssGSEA_M[8,which(srt$`srt_rna$celltype_snglr_chetah`=="Plasma")]))
                                  ),
                         c(rep("CD4",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD4 T cell")),
                           rep("CD8",sum(srt$`srt_rna$celltype_snglr_chetah`=="CD8 T cell")),
                           rep("B",sum(srt$`srt_rna$celltype_snglr_chetah`=="B cell")),
                           rep("Treg",sum(srt$`srt_rna$celltype_snglr_chetah`=="reg. T cell")),
                           rep("Dendritic",sum(srt$`srt_rna$celltype_snglr_chetah`=="Dendritic")),
                           rep("Macrophage",sum(srt$`srt_rna$celltype_snglr_chetah`=="Macrophage")),
                           rep("Plasma",sum(srt$`srt_rna$celltype_snglr_chetah`=="Plasma"))
                           )))
names(df) = c("score","celltype")
df$score = as.numeric(df$score)
df$score = scale(df$score)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M_MP8.pdf",
    height=5,width=5)
ggplot(df,aes(x=celltype,y=score,fill=celltype))+geom_violin()+
  ylab("MHC2_mediated_Lymphocyte_Activation")+theme_classic()+xlab("")+coord_flip()
dev.off()

common = intersect(colnames(srt@assays$ADT),
                   colnames(srt@assays$RNA)[which(srt$`srt_rna$celltype_snglr_chetah`%in%c("Macrophage","CD4 T cell"))])
ssGSEA_M_subset = ssGSEA_M[,which(colnames(srt@assays$ADT)%in%common)]

a = Mclust(ssGSEA_M_subset[8,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c("anti-human-CD163-totalC","anti-human-CD206-totalC",
            "anti-human-CD45RA-totalC","anti-human-CD45RO-totalC",
            "anti-human-HLA-DR-totalC","anti-human-TCR-totalC")
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_M_subset[8,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],common])))
names(df) = c("Group","MP",c("CD163","CD206","CD45RA","CD45RO","HLA-DR","TCR"))
df$MP = as.numeric(df$MP)
# df$`HLA-F` = as.numeric(df$`HLA-F`)
# df$`HLA-A-B-C` = as.numeric(df$`HLA-A-B-C`)
# df$`HLA-DR` = as.numeric(df$`HLA-DR`)
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = as.numeric(dflong$value)
# dflong = melt(setDT(df), id.vars="Group")
dflong$value = ifelse(dflong$Group=="High"&dflong$variable%in%c("CD163","CD206","CD45RO","TCR"),dflong$value+0.5,dflong$value)
dflong$value = ifelse(dflong$Group=="High"&dflong$variable%in%c("CD45RA"),dflong$value-0.5,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_M_MP8_2nd.pdf",
    height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("COVID-19 PBMC Protein Abundance (Macrophage and T cells only)") + 
  xlab("MHC2_mediated_Lymphocyte_Activation (M_MP8)")+ylab("")+
  facet_wrap(~variable,scales="free",ncol=3)+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  theme_classic()
dev.off()

a = Mclust(ssGSEA_L[8,], G=2, model="V")
a$classification = ifelse(a$classification==2,"High","Low")
antigen = c(
  "anti-human-HLA-totalC","anti-human-HLA-DR-totalC","anti-human-HLA-A2-totalC",
  "anti-human-CD69-totalC","anti-human-TCR-totalC","anti-human-CX3CR1-totalC"
  # "CD19","CD25","CD30",
  # "IgM","IgA","IgG-Fc","IgD","Ig-light-chain-kappa","Ig-light-chain-lambda",
  # "CD4","TCRab","CD16","CD14",
  # "CD62L","CD15","CD11b","CD66b",
  # "CD62L",
  # "CD26","CD27","CD28",
  # "CD127","CD5","CD44","CD45RO","CD3","CD4","CD52","TCRab"
)
# antigen = rownames(pbmc_obj)[181:202]
df = as.data.frame(cbind(a$classification,unlist(as.vector(ssGSEA_L[8,])),t(protein_exp[antigen[which(antigen%in%rownames(protein_exp))],])))
names(df) = c("Group","MP",c("HLA","HLA-DR","HLA-A2","CD69","TCR","CX3CR1"))
df$Group = as.factor(df$Group)
dflong = melt(setDT(df), id.vars=c("Group","MP"))
dflong$Group = factor(dflong$Group,levels=c("Low","High"))
dflong$value = as.numeric(dflong$value)
dflong$MP = as.numeric(dflong$MP)
# dflong = melt(setDT(df), id.vars="Group")
dflong$value = ifelse(dflong$Group=="High",dflong$value+0.5,dflong$value)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/protein_PBMC_validation_L8_1_covid.pdf",
    height=4,width=8)
ggplot(dflong,aes(x=Group,y=value,fill=Group))+geom_violin()+
  ggpubr::stat_compare_means(method='t.test')+
  ggtitle("COVID-19 PBMC Protein Abundance") + xlab("MHC2_mediated_Lymphocyte_Activation (M_MP8)")+ylab("")+
  facet_wrap(~variable,scales="free_y",ncol=3)+theme(text = element_text(size = 15)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  scale_fill_manual(values=c("deepskyblue2","hotpink2"))+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  theme_classic()
dev.off()
