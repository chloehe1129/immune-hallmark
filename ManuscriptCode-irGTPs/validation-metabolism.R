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
library(ggplot2)
library(ggpubr)
# library(igraph)
# library(ComplexHeatmap)
# library(org.Hs.eg.db)
# library(msigdbr)
# library(fgsea)
# library(clusterProfiler)
library(Matrix)
library(reshape2)
library(data.table)
library(AUCell)

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
gs_list = c(gs_list_L,gs_list_M)

######################### T metabolism CD8+ activation stages ######################### 

mat = readMM("/Volumes/she4/hallmark/data/validation-data/metabolism_T/GSE211602_matrix_AllSamples.mtx")
gene = read.table("/Volumes/she4/hallmark/data/validation-data/metabolism_T/GSE211602_genes_AllSamples.tsv")
cell = read.table("/Volumes/she4/hallmark/data/validation-data/metabolism_T/GSE211602_barcodes_AllSamples.tsv")
colnames(mat) = cell$V1
rownames(mat) = toupper(gene$V2)
meta = read.table("/Volumes/she4/hallmark/data/validation-data/metabolism_T/GSE211602_CellMetaData_AllSamples.txt",header = T)
rownames(meta) = meta$CellBarcode
srt = CreateSeuratObject(mat,meta.data = meta)
srt$CellState = ifelse(srt$CellState=="Na\xefve","Naive",srt$CellState)

gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
H_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(H_list) = unique(gene_sets$gs_name)

cells_AUC <- AUCell_run(srt@assays$RNA$counts, geneSets=list(geneset = gs_list_L[[9]]))
cells_AUC_gly <- AUCell_run(srt@assays$RNA$counts, geneSets=list(geneset = H_list[[19]]))
cells_AUC_ox <- AUCell_run(srt@assays$RNA$counts, geneSets=list(geneset = H_list[[36]]))
srt$metabolism_score = cells_AUC@assays@data@listData[["AUC"]]
srt$glycolysis = cells_AUC_gly@assays@data@listData[["AUC"]]
srt$oxphos = cells_AUC_ox@assays@data@listData[["AUC"]]
df = srt@meta.data[,c("UMAPx","UMAPy","metabolism_score","glycolysis","oxphos","CellState")]
df = df[which(df$CellState%in%c("Naive","Early","Mid","Mid-Late","Late")),]
p1=ggplot(df,aes(x=UMAPx,y=UMAPy,col=CellState))+geom_point()+theme_classic()
p2=ggplot(df,aes(x=UMAPx,y=UMAPy,col=metabolism_score))+geom_point()+
  scale_color_gradient2(midpoint = 0.15, low = "blue", mid = "white",
                          high = "red", space = "Lab" )+ggtitle(names(gs_list_L)[9])+theme_classic()+ 
  guides(col = guide_colourbar(title = "GS Density"))
p3=ggplot(df,aes(x=UMAPx,y=UMAPy,col=glycolysis))+geom_point()+
  scale_color_gradient2(midpoint = 0.04, low = "blue", mid = "white",
                        high = "red", space = "Lab" )+ggtitle("Glycolysis")+theme_classic() +
guides(col = guide_colourbar(title = "GS Density"))
p4=ggplot(df,aes(x=UMAPx,y=UMAPy,col=oxphos))+geom_point()+
  scale_color_gradient2(midpoint = 0.15, low = "blue", mid = "white",
                        high = "red", space = "Lab" )+ggtitle("Oxidative Phosphorylation")+theme_classic() +
guides(col = guide_colourbar(title = "GS Density"))

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Metabolism_T_state.pdf",
    height=3.5,width=15)
print(patchwork::wrap_plots(p1,p2,p3,p4,ncol=4))
dev.off()

# pdf("~/desktop/T_state.pdf",onefile = T,height=7,width=14)
# for(i in 1:19){
#   cells_AUC <- AUCell_run(srt@assays$RNA$counts, geneSets=list(geneset = gs_list_L[[i]]))
#   srt$metabolism_score = cells_AUC@assays@data@listData[["AUC"]]
#   df = srt@meta.data[,c("UMAPx","UMAPy","metabolism_score","CellState")]
#   df = df[which(df$CellState%in%c("Naive","Early","Mid","Mid-Late","Late")),]
#   p1=ggplot(df,aes(x=UMAPx,y=UMAPy,col=CellState))+geom_point()
#   p2=ggplot(df,aes(x=UMAPx,y=UMAPy,col=metabolism_score))+geom_point()+
#     scale_color_gradient2(midpoint = 0.15, low = "blue", mid = "white",
#                           high = "red", space = "Lab" )+ggtitle(names(gs_list_L)[i])
#   print(p1+p2)
# }
# dev.off()


######################### NK metabolism lili data ######################### 

nk<-readRDS(paste0("/Volumes/she4/hallmark-sc/data/","validation-data/lili/nk.rds"))
cells_AUC_meta <- AUCell_run(nk@assays$RNA$counts, geneSets=list(Metabolism = gs_list_L[[9]]))
nk$metabolism_score = cells_AUC_meta@assays@data@listData[["AUC"]]
cells_AUC_gly <- AUCell_run(nk@assays$RNA$counts, geneSets=list(glycolysis = H_list[[19]]))
cells_AUC_ox <- AUCell_run(nk@assays$RNA$counts, geneSets=list(oxphos = H_list[[36]]))
nk$glycolysis = cells_AUC_gly@assays@data@listData[["AUC"]]
nk$oxphos = cells_AUC_ox@assays@data@listData[["AUC"]]
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Metabolism_NK_state.pdf",
    height=3.5,width=12)
FeaturePlot(nk,features = c("metabolism_score","glycolysis","oxphos"),cols = c("cornflowerblue", "hotpink2"),
            ncol=3)
dev.off()


########################################################################################
# # Antiviral Response
# srt <- readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu.rds"))
# srt = UpdateSeuratObject(srt)
# srt = subset(srt,downsample=100)
# 
# srt@meta.data$drug = ifelse(srt@meta.data$Donor%in%c("C1","C3","C4"),
#                             "Azithromycin",ifelse(srt@meta.data$Donor%in%c("H1","H2","H3","H4","H5","H6"),"Healthy","Remdesivir"))
# 
# ssGSEA = corto::ssgsea(inmat=as.matrix(srt@assays$RNA@counts),groups=gs_list,minsize = 0)
# ssGSEA = as.data.frame(t(ssGSEA))
# dat = as.data.frame(cbind(ssGSEA[,1],ssGSEA[,22],ssGSEA[,23],srt$drug))
# names(dat) = c("Antiviral_Defense_Network","Interferon_induced_Antiviral_Defense","Type.I_Interferon","Intervention")
# dat = melt(setDT(dat),id.vars = "Intervention")
# 
# print("done")
# 
# dat$value = as.numeric(dat$value)
# pdf(paste0(path_result,"validation/antiviral_COVID_ssgsea_scale.pdf"),width=8,height=3)
# ggplot(data=dat[which(dat$Intervention!="Azithromycin"),],
#        aes(x=Intervention,y=value,fill=Intervention))+geom_violin()+
#   facet_wrap(~variable,scale='free_y',ncol=3)+stat_compare_means(method="t.test",label="p.signif")+
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
# dev.off()
# 
# # score_ssGSEA = NULL
# # for (i in 1:length(unique(srt$drug))){
# #   score_subset = ssGSEA[which(srt$drug==unique(srt$drug)[i]),]
# #   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean))
# # }
# # rownames(score_ssGSEA) = unique(srt$drug)
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # score_ssGSEA = range01(score_ssGSEA)
# # dat = as.data.frame(cbind(unname(score_ssGSEA[,1]),unique(srt$drug),rep("Antiviral Defense",3)))
# # dat = dat[-c(3,6),]
# # dat$V1 = round(as.numeric(dat$V1),3)
# # names(dat)[3] = "Pathway"
# # names(dat)[2] = "Drug"
# # pdf(paste0(path_result,"validation/antiviral_COVID_ssgsea_scale.pdf"),width=15)
# # Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "COVID-19",
# #         column_names_max_height = unit(12, "cm"),
# #         cell_fun = function(j, i, x, y, width, height, fill) {
# #           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
# #         })
# # dev.off()
# 
# # ce <- compute.mca(object = srt)
# # res <- compute.kld(coembed = ce,
# #                    genes.use = intersect(rownames(ce),rownames(srt)), 
# #                    n.grids = 100,
# #                    gene.set.list = gs_list_M,
# #                    gene.set.cutoff = 3,
# #                    n.times = 100)
# # cells <- colnames(srt)
# # 
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# # 
# # res = matrix(0,ncol(srt),9)
# # for (gs in 1:9){
# #   res[,gs] <- run.rwr(el = el, gene_set = gs_list_M[[gs]], cells = cells)
# # }
# # write.csv(res,paste0(path_result,"validation/antiviral_M_GSD_scale.csv"))
# 
# 
# # pdf(paste0(path_result,"validation/antiviral_covid.pdf"),height = 10, width = 22)
# # DimPlot(srt,
# #         group.by = c("Admission","Ventilated","Status","cell.type.fine","drug"),
# #         label=T,
# #         raster = T)
# # dev.off()
# # 
# # Idents(srt) = "cell.type.fine"
# # 
# # # Plot ISG module score
# # ISG = read.delim(paste0(path_data,"validation-data/antiviral/covid/ISG_genes.txt"),header = F)
# # ISG = list(ISG$V1[28:47])
# # names(ISG) = "ISG"
# # ce <- compute.mca(object = srt)
# # cells <- colnames(srt)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# # cv <- run.rwr(el = el, gene_set = ISG[[1]], cells = cells)
# # srt@meta.data$ISG_module_score <- cv[colnames(srt)]
# # pdf(paste0(path_result,"validation/antiviral_GSDensity_ISG_covid.pdf"),height = 7, width = 7)
# # FeaturePlot(srt,features = "ISG_module_score",label=T,raster = T)
# # dev.off()
# # 
# # # Plot exhaustion module score
# # # exh = read.delim(paste0(path_data,"validation-data/antiviral/covid/exhaustion.txt"),header = T)
# # exh = c("TIGIT","TOX","PDCD1","CD160","CD244","CTLA4","BTLA","HAVCR2","LAG3")
# # exh = list(exh)
# # names(exh) = "Exhaustion"
# # ce <- compute.mca(object = srt)
# # cells <- colnames(srt)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# # cv <- run.rwr(el = el, gene_set = exh[[1]], cells = cells)
# # srt@meta.data$Exhaustion_module_score <- cv[colnames(srt)]
# # pdf(paste0(path_result,"validation/antiviral_GSDensity_Exhaustion_covid.pdf"),height = 7, width = 7)
# # FeaturePlot(srt,features = "Exhaustion_module_score",label=T,raster = T)
# # dev.off()
# # 
# # annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"),row.names = 1)
# # gs_list_L = list()
# # file_names <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
# # for (i in 1:length(file_names)){
# #   gs_list_L[[i]] = read.csv(file_names[i],row.names = 1)[,1]
# # }
# # names(gs_list_L) = names(annotation)
# # 
# # annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"),row.names = 1)
# # gs_list_M = list()
# # file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
# # for (i in 1:length(file_names)){
# #   gs_list_M[[i]] = read.csv(file_names[i],row.names = 1)[,1]
# # }
# # names(gs_list_M) = names(annotation)
# 
# # ce <- compute.mca(object = srt)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # cells <- colnames(srt)
# # el <- compute.nn.edges(coembed = ce, nn.use = 600)
# # 
# # for (gs in 1:13){
# #   cv <- run.rwr(el = el, gene_set = gs_list_M[[gs]], cells = cells)
# #   compute.cell.label <- function(cell_vec){
# #     nm <- names(cell_vec)
# #     m <- locmodes(cell_vec, mod0 = 4)
# #     split.value <- m$locations[3]
# #     cell.label <- ifelse(cell_vec < split.value, "negative", "positive")
# #     return(cell.label)
# #   }
# #   #cl = ifelse(cv<mean(cv),"negative","positive")
# #   #cl <- compute.cell.label(cv)
# #   srt@meta.data$geneset <- cv[colnames(srt)]
# #   #srt@meta.data$geneset_bin <- cl[colnames(srt)]
# #   colnames(srt@meta.data)[which(colnames(srt@meta.data)=="geneset")] = names(gs_list_M[gs])
# #   
# #   assign(paste0("p",gs),FeaturePlot(srt,
# #                                 features = names(gs_list_M[gs]),
# #                                 label=T,
# #                                 raster = T))
# # }
# # pdf(paste0(path_result,"validation/antiviral_GSDensity_M_covid.pdf"),height = 20, width = 25)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,
# #              nrow=4,ncol=4)
# # dev.off()
# # 
# # ce <- compute.mca(object = srt)
# # res <- compute.kld(coembed = ce,
# #                    genes.use = intersect(rownames(ce),rownames(srt)), # this intersection is to select only genes, not cells.
# #                    n.grids = 100,
# #                    gene.set.list = gs_list_L,
# #                    gene.set.cutoff = 3,
# #                    n.times = 100)
# # cells <- colnames(srt)
# # el <- compute.nn.edges(coembed = ce, nn.use = 600)
# # 
# # for (gs in 1:21){
# #   cv <- run.rwr(el = el, gene_set = gs_list_L[[gs]], cells = cells)
# #   compute.cell.label <- function(cell_vec){
# #     nm <- names(cell_vec)
# #     m <- locmodes(cell_vec, mod0 = 4)
# #     split.value <- m$locations[3]
# #     cell.label <- ifelse(cell_vec < split.value, "negative", "positive")
# #     return(cell.label)
# #   }
# #   #cl = ifelse(cv<mean(cv),"negative","positive")
# #   #cl <- compute.cell.label(cv)
# #   srt@meta.data$geneset <- cv[colnames(srt)]
# #   colnames(srt@meta.data)[which(colnames(srt@meta.data)=="geneset")] = names(gs_list_L[gs])
# #   #srt@meta.data$geneset_bin <- cl[colnames(srt)]
# # 
# #   assign(paste0("p",gs),FeaturePlot(srt,
# #                                     features = names(gs_list_L[gs]),
# #                                     label=T,
# #                                     raster = T))
# # }
# # pdf(paste0(path_result,"validation/antiviral_GSDensity_L_covid.pdf"),height = 15, width = 35)
# # grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,
# #              nrow=3,ncol=7)
# # dev.off()
# 
# # Antiviral Defense through interacting with viral RNA
# 
# # srt_subset = subset(srt,subset=drug%in%c("Azithromycin","Remdesivir"))
# # ssGSEA = corto::ssgsea(as.matrix(srt_subset@assays$RNA@counts),groups=gs_list_M)
# # ssGSEA = as.data.frame(t(ssGSEA))
# # score_ssGSEA = NULL
# # for (i in 1:length(unique(srt_subset$drug))){
# #   score_subset = ssGSEA[which(srt_subset$drug==unique(srt_subset$drug)[i]),]
# #   score_ssGSEA = rbind(score_ssGSEA,apply(score_subset,2,mean))
# # }
# # rownames(score_ssGSEA) = unique(srt_subset$drug)
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # score_ssGSEA = range01(score_ssGSEA)
# # pdf(paste0(path_result,"validation/antiviral_ssGSEA.pdf"),height=7,width=10)
# # Heatmap(as.matrix(score_ssGSEA),name="MP ssGSEA score", cluster_columns = F,column_title = "COVID 19 Atlas",
# #         column_names_max_height = unit(12, "cm"),
# #         cell_fun = function(j, i, x, y, width, height, fill) {
# #           grid.text(sprintf("%.2f", score_ssGSEA[i, j]), x, y, gp = gpar(fontsize = 10))
# #         })
# # dev.off()
# 
# # ce <- compute.mca(object = srt_subset)
# # cells <- colnames(srt_subset)
# # el <- compute.nn.edges(coembed = ce, nn.use = 1000)
# # source("/rsrch4/home/bcb/she4/hallmark-sc/code/GSDensity_code.R")
# # cv <- run.rwr(el = el, gene_set = gs_list_M[[1]], cells = cells)
# # cl = ifelse(cv<mean(cv),"Azithromycin","Remdesivir")
# # cm = confusionMatrix(as.factor(srt_subset@meta.data$drug),as.factor(cl))
# # print(cm)
# # azi_score = mean(cv[which(srt_subset$drug=="Azithromycin")])
# # rem_score = mean(cv[which(srt_subset$drug=="Remdesivir")])
