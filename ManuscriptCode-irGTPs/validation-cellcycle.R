path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

library(Seurat)
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

srt = readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD-L.rds"))
pdf(paste0(path_result,"validation/CellCycle_Lymphocyte_all_COVID19_UMAP.pdf"),height = 26, width = 32)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,nrow=4,ncol=5)
dev.off()

p1=VlnPlot(srt,features="Cell_Cycle_Progression",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  stat_compare_means(size=8, label.y = 0.00015,label.x = 4)+xlab("")+ggtitle("L_MP1: Cell Cycle 1")
p5=VlnPlot(srt,features="CellCycle_Mitosis",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+theme(text = element_text(size = 30))+ 
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  stat_compare_means(size=8, label.y = 0.00013,label.x = 4)+xlab("")+ggtitle("L_MP6: Cell Cycle 2")
p2=VlnPlot(srt,features="CellCycle_Immune_Response",pt.size = 0,
           cols=c("slateblue","hotpink2","brown4","plum1","palegreen4","cornflowerblue","goldenrod1","darkorange3","skyblue"))+theme(text = element_text(size = 30))+ 
  theme(text = element_text(size = 30))+ 
  theme(axis.text.x = element_text(size=25,angle = 90, hjust = 1))+
  theme(axis.text.y = element_text(size=25,angle = 0, hjust = 1))+
  stat_compare_means(size=8, label.y = 0.00013,label.x = 4)+xlab("")+ggtitle("L_MP11: Cell Cycle Immune Response")

pdf(paste0(path_result,"validation/CellCycles_COVID19_Violin.pdf"),height = 15,width=14)
grid.arrange(p1,p5,p2,nrow=3,ncol=1)
dev.off()
pdf(paste0(path_result,"validation/CellCycles_COVID19_Violin_horizontal.pdf"),
    height = 6,width=30)
grid.arrange(p1,p5,p2,nrow=1,ncol=3)
dev.off()

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

files = list.files("/Volumes/she4/hallmark/data/validation-data/cellcycle",full.names = T)
options("Seurat.object.assay.version" = "v3")
mat1 = read.table(files[1])
srt1 = CreateSeuratObject(mat1)
for(i in 2:4){
  mat1 = read.table(files[i])
  seu = CreateSeuratObject(mat1)
  srt1 = merge(srt1,seu)
  print(i)
}
srt1<- JoinLayers(srt1,  assay = "RNA")
srt1 = NormalizeData(srt1)
srt1 = ScaleData(srt1)
srt1 = FindVariableFeatures(srt1)
srt1 = RunPCA(srt1)
srt1 = FindNeighbors(srt1)
srt1 = FindClusters(srt1)
srt1 = RunUMAP(srt1, dims = 1:10, verbose = F)

# srt = readRDS(paste0(path_data,"validation-data/antiviral/covid/blish_covid.seu-GSD-L.rds"))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
options("Seurat.object.assay.version" = "v3")
srt1 <- CellCycleScoring(srt1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

ssGSEA = corto::ssgsea(inmat=as.matrix(srt1@assays$RNA$counts),groups=gs_list_L)
ssGSEA = as.data.frame(t(ssGSEA))
rownames(ssGSEA) = colnames(srt1)
srt1 = AddMetaData(srt1,metadata=ssGSEA)

p1=VlnPlot(srt1,features="Cell_Cycle_1",group.by = "Phase",pt.size = 0,cols=c("skyblue","hotpink2","goldenrod1"))+ 
  stat_compare_means(comparisons=list(c("G2M", "S"),c("G1", "S"),c("G1", "G2M")),label = "p.signif")+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylim(-2, 7)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+geom_boxplot()

p2=VlnPlot(srt1,features="Cell_Cycle_2",group.by = "Phase",pt.size = 0,cols=c("skyblue","hotpink2","goldenrod1")) + 
  stat_compare_means(comparisons=list(c("G2M", "S"),c("G1", "S"),c("G1", "G2M")),label = "p.signif")+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylim(-2, 10)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+geom_boxplot()

p3=VlnPlot(srt1,features="CellCycle_Immune_Response",group.by = "Phase",pt.size = 0,cols=c("skyblue","hotpink2","goldenrod1")) + 
  stat_compare_means(comparisons=list(c("G2M", "S"),c("G1", "S"),c("G1", "G2M")),label = "p.signif")+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ylim(-2, 6)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+geom_boxplot()
pdf("/Volumes/she4/hallmark/results/validation/CellCycles_3stages.pdf",height = 12,width=7)
patchwork::wrap_plots(p1,p2,p3,nrow=3)
dev.off()
