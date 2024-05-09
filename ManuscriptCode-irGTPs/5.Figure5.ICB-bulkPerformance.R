library(corto)
library(umap)
library(ggplot2)
library(ggforce)
library(Matrix)

path_code = "/Volumes/she4/hallmark/code/"
path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_code = "/rsrch4/home/bcb/she4/hallmark/code/"
path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

# Our MPs
annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
gs_list_L = list()
file_names_L <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names_L)){
  gs_list_L[[i]] = read.csv(file_names_L[i],row.names = 1)[,1]
}
names(gs_list_L) = paste0("L_",names(annotation))

annotation = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))
gs_list_M = list()
file_names_M <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
# file_names <- dir(paste0(path_result,"cluster_methodTirosh2"),full.names = TRUE,pattern = "_based") #where you have your files
for (i in 1:length(file_names_M)){
  gs_list_M[[i]] = read.csv(file_names_M[i],row.names = 1)[,1]
}
names(gs_list_M) = paste0("M_",names(annotation))

gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
H_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(H_list) = unique(gene_sets$gs_name)

kegg_sets = read.csv(paste0(path_data,"kegg_immune_sets.csv"))[,-1]
KEGG_list = list()
for (i in 1:length(unique(kegg_sets$name))){
  KEGG_list[[i]] = unique(kegg_sets[which(kegg_sets$name==unique(kegg_sets$name)[i]),"symbol"])
}
names(KEGG_list) = unique(kegg_sets$name)

gs_list = c(gs_list_L,gs_list_M,H_list,KEGG_list)

############## LUNG ICB ############## 
mat = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE111414_lung_ICB_RNA/GSE111414_gene_counts.csv")
mat = mat[!duplicated(mat$SYMBOL),]
mat = mat[!is.na(mat$SYMBOL),]
rownames(mat) = mat$SYMBOL
mat = mat[,-c(1:4)]
mat = mat[,names(mat)[grep("N1",names(mat))]]
ssgsea = corto::ssgsea(mat,gs_list,minsize = 0)
colnames(ssgsea) = colnames(mat)
ssgsea = as.data.frame(t(ssgsea))
meta = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE111414_lung_ICB_RNA/meta.csv")
ssgsea = ssgsea[meta$patient,]

############## Bladder ICB ############## 
mat = read.table("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE176307_uro/GSE176307_baci_rsem_RS_BACI_headers_tab.txt",
                 sep="\t",row.names = 1,header = T)
mat = mat[-1,]
genes = rownames(mat)
mat = apply(mat,2,as.numeric)
rownames(mat) = genes
ssgsea = corto::ssgsea(mat,gs_list,minsize = 0)
colnames(ssgsea) = colnames(mat)
ssgsea = as.data.frame(t(ssgsea))
meta = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE111414_lung_ICB_RNA/meta.csv")
ssgsea = ssgsea[meta$patient,]

############## BREAST ICB ############## 
# mat = readMM("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE169246_breast_ICB_scRNA/GSE169246_TNBC_RNA.counts.mtx.gz")
# genes = read.table("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE169246_breast_ICB_scRNA/TNBC_RNA.feature.tsv")
# barcode = read.table("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE169246_breast_ICB_scRNA/TNBC_RNA.barcode.tsv")
# barcode = as.data.frame(barcode[grep("Pre",barcode$V1),])
# barcode = as.data.frame(barcode[grep("_b",barcode$`barcode[grep("Pre", barcode$V1), ]`),])
# names(barcode) = "V1"
# colnames(mat) = barcode$V1
# rownames(mat) = genes$V1
# srt_sc = CreateSeuratObject(mat)
# meta = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE168204_melanoma_ICB_RNA/meta_melanoma.csv")
# srt_sc$response = NA
# for(i in 1:nrow(srt_sc)){
#   srt_sc$response[i] = meta[which(meta$Sample == rownames(srt_sc)[i]),"Resp_NoResp"]
# }
# srt_bulk = AggregateExpression(srt_sc,assays ="RNA",
#                                group.by='response',
#                                return.seurat = T)
# ssgsea = corto::ssgsea(as.matrix(srt_bulk@assays$RNA$scale.data),gs_list,minsize = 0)
# ssgsea = as.data.frame(t(ssgsea))

############## Melanoma ICB ############## 
mat = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE168204_melanoma_ICB_RNA/GSE168204_MGH_counts.csv")
mat = mat[!duplicated(mat$Gene),]
rownames(mat) = mat$Gene
mat = mat[,-1]
ssgsea = corto::ssgsea(mat,gs_list,minsize = 0)
colnames(ssgsea) = colnames(mat)
ssgsea = as.data.frame(t(ssgsea))
rownames(ssgsea) = gsub("\\.","_",rownames(ssgsea))
meta = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE168204_melanoma_ICB_RNA/meta_melanoma.csv")
ssgsea = ssgsea[meta$Sample,]

############## Liver ICB ############## 
mat = read.table("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE140901_liver/GSE140901_processed_data.txt")
colnames(mat) = mat[1,]
mat = mat[-1,]
rownames(mat) = mat$Gene
mat = mat[,-1]
genes = rownames(mat)
mat = apply(mat,2,as.numeric)
rownames(mat) = genes 
ssgsea = corto::ssgsea(mat,gs_list,minsize = 0)
colnames(ssgsea) = colnames(mat)
ssgsea = as.data.frame(t(ssgsea))
meta = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/data/figure5_validation/GSE140901_liver/meta.csv")
ssgsea = ssgsea[meta$patient,]

############## ssgsea and check separation ############## 
ssgsea_MP = cbind(ssgsea[,c("M_Eosinophil_chemotaxis",
                                        "L_CellCycle_Immune_Response",
                                        "L_Cell_Cycle_1",
                                        "L_Cell_Adhesion_2.CD52..",
                                        "M_MHC2_mediated_Lymphocyte_Activation",
                                        "L_Lipid_Localization_TCR_synapse",
                                        "L_Cytokine_Receptor_Signaling",
                                        "L_Antigen_Presentation",
                                        "L_MHC.II_mediated_immunity",
                                        "L_Cell_Adhesion_1")])
# iter = sample(3686:4444,1, replace=F)  #2732(melanoma), 3686(liver)
set.seed(88)
ssgsea_MP[is.na(ssgsea_MP)]=0
ss = umap(as.matrix(ssgsea_MP),n_neighbors=10)
km.final <- kmeans(ssgsea_MP, 2)
def = as.data.frame(cbind(ss$layout,meta$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p1=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE111414 - 10 irMPs")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")
p1

ssgsea_H = cbind(ssgsea[,c(29:78)])
set.seed(88)
ssgsea_H[is.na(ssgsea_H)]=0
ss = umap(ssgsea_H,n_neighbors=10)
km.final <- kmeans(ssgsea_H, 2)
def = as.data.frame(cbind(ss$layout,meta$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p2=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE111414 - 50 Hallmarks")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")

ssgsea_K = cbind(ssgsea[,c(79:107)])
set.seed(88)
ssgsea_K[is.na(ssgsea_K)]=0
ss = umap(ssgsea_K,n_neighbors=10)
km.final <- kmeans(ssgsea_K, 2)
def = as.data.frame(cbind(ss$layout,meta$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p3=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE111414 - 29 irKEGGs")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.GSE111414_UMAP.pdf",
    height=4,width=12)
p1 + p2 + p3 + patchwork::plot_layout(guides = "collect",ncol=3) & 
  theme(legend.position = 'bottom',legend.justification = "left")
dev.off()