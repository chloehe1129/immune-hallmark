# library(stringi)
library(Seurat)
library(dplyr)
# library(CellChat)
# library(AUCell)

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

################ GSE120575 ################

cell_annot = read.csv(paste0(path_data,"/validation-data/GSE120575/cell.annotations.csv"),row.names = 1)
CD8_cell_cluster = read.csv(paste0(path_data,"/validation-data/GSE120575/CD8cluster_cell_name_mapping.csv"),row.names = 1)
CD8_cell_cluster$CD8_subtype = ifelse(CD8_cell_cluster$Cluster=="CD8_1","Exhaustion/CellCycle",
                                      ifelse(CD8_cell_cluster$Cluster=="CD8_2","Exhaustion/HSP",
                                             ifelse(CD8_cell_cluster$Cluster=="CD8_3","Exhaustion",
                                                    ifelse(CD8_cell_cluster$Cluster%in%c("CD8_4","CD8_6"),"Memory/Effector",
                                                           "Early.Activated.Cells"))))
cell_annot$CD8_subtype = ""
for(i in 1:nrow(cell_annot)){
  if(length(which(rownames(CD8_cell_cluster)==cell_annot$title[i]))!=0){
    cell_annot$CD8_subtype[i] = CD8_cell_cluster[which(rownames(CD8_cell_cluster)==cell_annot$title[i]),"CD8_subtype"]
  }
}

cell_annot$timepoint = unlist(stri_split_fixed(str = cell_annot$characteristics..patinet.ID..Pre.baseline..Post..on.treatment., pattern = "_", n = 2))[seq(1,nrow(cell_annot)*2,2)]
rownames(cell_annot) = cell_annot$title
mat = read.csv(paste0(path_data,"validation-data/GSE120575/counts.csv"),row.names = 1,header=T)
colnames(mat) = cell_annot$title
mat[is.na(mat)] = 0
print(sum(is.na(mat)))

srt = CreateSeuratObject(mat,meta.data = cell_annot)
saveRDS(srt,paste0(path_data,"/validation-data/GSE120575/seurat.rds"))

srt = readRDS("/Volumes/she4/hallmark/data/validation-data/GSE120575/seurat.rds")
doublet = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE120575/doublet_score.csv",row.names = 1)
doublet$doublet = ifelse(doublet$X0>unname(quantile(doublet$X0, probs = 0.95)),"doublet","singlet")
rownames(doublet) = colnames(srt)
srt = AddMetaData(srt,metadata =  doublet)

srt = NormalizeData(srt)
srt = FindVariableFeatures(srt)
srt = ScaleData(srt)
srt = RunPCA(srt)
srt = FindNeighbors(srt)
srt = FindClusters(srt)
srt = RunUMAP(srt,dims=c(1:10))
Idents(srt) = "doublet"
saveRDS(srt,paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet.rds"))

###### (1) Re-annotate - Cell type proportion before and after ######

srt = readRDS(paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet.rds"))
celltypist = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE120575/celltypist.csv",row.names = 1)
srt = AddMetaData(srt,metadata = celltypist)
Idents(srt) = "majority_voting"
srt$reannotation = ifelse(srt$majority_voting%in%c("CRTAM+ gamma-delta T cells","gamma-delta T cells"),"gd T cells",
                          ifelse(srt$majority_voting%in%c("Intermediate macrophages","Macrophages"),"Macrophages",
                                 ifelse(srt$majority_voting%in%c("Non-classical monocytes","Classical monocytes"),"Monocytes",
                                 srt$majority_voting)))
DimPlot(srt,group.by = "reannotation",label=T,repel = T)
Idents(srt) = "reannotation"
saveRDS(srt,paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet_celltypist.rds"))

############## cell type composition before and after comparison ##############

df = as.data.frame(srt@meta.data[which(srt@meta.data$timepoint=="Post"),c("orig.ident","timepoint","cell.types","characteristics..response")])
df_pre = df[which(df$characteristics..response=="Responder"),]
res1 = NULL
for (i in 1:length(unique(df_pre$cell.types))){
  res1 = rbind(res1,c(length(which(df_pre$cell.types==unique(df_pre$cell.types)[i])),nrow(df_pre)))
}
res1 = as.data.frame(res1)
res1$`Cell Type` = unique(df_pre$cell.types)
res1$response = "Responder"

df_post = df[which(df$characteristics..response=="Non-responder"),]
res2 = NULL
for (i in 1:length(unique(df_post$cell.types))){
  res2 = rbind(res2,c(length(which(df_post$cell.types==unique(df_post$cell.types)[i])),nrow(df_post)))
}
res2 = as.data.frame(res2)
res2$`Cell Type` = unique(df_post$cell.types)
res2$response = "Non-responder"

res = as.data.frame(rbind(res1,res2))
names(res) = c("Abundance","Total","Cell Type","Response")
res = res[order(res$`Cell Type`),]

pval1 = NULL
for (ct in 1:length(unique(res$`Cell Type`))){
  dat <- data.frame(
    "celltype" = res[which(res$`Cell Type`==unique(res$`Cell Type`)[ct]),1],
    "non_celltype" = c(sum(res[which(res$`Cell Type`!=unique(res$`Cell Type`)[ct]&res$Response=="Responder"),1]),
                       sum(res[which(res$`Cell Type`!=unique(res$`Cell Type`)[ct]&res$Response=="Non-responder"),1])),
    row.names = c("Pre", "After"),
    stringsAsFactors = FALSE
  )
  pval1 = c(pval1,fisher.test(dat)$p.value)
}
names(pval1) = unique(res$`Cell Type`)
pval1 = p.adjust(pval1, method = "fdr", n = length(pval1))

res$Response <- factor(res$Response, levels = c("Responder", "Non-responder"))
stat.test = as.data.frame(cbind(names(pval1),rep("Relative Abundance",11),rep("Pre",11),rep("Post",11),unname(pval1)))
names(stat.test) = c("Cell Type",".y.","group1","group2","p.adj")
stat.test = as_tibble(stat.test)
stat.test$`Cell Type`=as.factor(stat.test$`Cell Type`)
stat.test$p.adj.sig = ifelse(as.numeric(stat.test$p.adj)<0.00001,"****",
                             ifelse(as.numeric(stat.test$p.adj)<0.0001,"***",
                                    ifelse(as.numeric(stat.test$p.adj)<0.001,"**",
                                           ifelse(as.numeric(stat.test$p.adj)<0.05,"*","ns"))))
res$`Relative Abundance` = res$Abundance/res$Total
bp <- ggbarplot(
  data = res, x = "Cell Type", y = "Relative Abundance",
  fill= "Response", palette = c("cornflowerblue", "pink2"),
  position = position_dodge(0.8)
)
# Add p-values onto the bar plots
stat.test <- stat.test %>%
  rstatix::add_xy_position(data=res,formula=`Relative Abundance`~Response, x = "dose", dodge = 0.8)
stat.test$xmin = (seq(1,22,2)+0.3)/2
stat.test$xmax = (seq(2,22,2)+0.3)/2
stat.test$p.adj = round(as.numeric(stat.test$p.adj),2)
p2 = bp + stat_pvalue_manual(hide.ns = FALSE,
                             stat.test,  label = "p.adj.sig",
                             tip.length = 0.01,y.position = 0.4,bracket.size = 0.6
)+scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
p2
pdf("/Volumes/she4/hallmark/results/validation/Figure5_scRNA_beforeVafter.pdf",height=12,width=12)
wrap_plots(p1,p2,ncol=1)
dev.off()

# markers = FindAllMarkers(srt,assay="RNA",only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
# markers = markers[which(markers$p_val_adj<0.05),]
# markers = as.data.frame(markers%>%
#                           group_by(cluster)%>%
#                           slice_max(n=50,order_by=avg_log2FC))
#
# srt_subset = subset(srt,subset = majority_voting%in%c("CD16- NK cells","CD16+ NK cells",
#                                                 "Regulatory T cells","Tcm/Naive helper T cells","Tem/Effector helper T cells",
#                                                 "Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells"))
# srt_subset_subset = subset(srt_subset,subset=timepoint=="Post")
# df = as.data.frame(t(table(srt_subset_subset$characteristics..response,srt_subset_subset$majority_voting)))
# df = reshape2::acast(df, Var1~Var2, value.var="Freq")
# prop = apply(df,2,prop.table)
# prop = melt(prop)
# names(prop) = c("CellType","Sample","Proportion")
# p1=ggplot(prop, aes(fill=CellType, y=Proportion, x=Sample)) +
#   geom_bar(position="fill", stat="identity")+
#   ggtitle("Cell Type Abundance for each sample")+
#   scale_fill_manual(values=c("skyblue","blue",
#                                       "mediumpurple1",
#                                       "limegreen","lightgreen",
#                                       "brown","pink"))
# wrap_plots(p1+p2,ncol=2)

# srt_subset_subset = subset(srt,subset=characteristics..response=="Non-responder")
# df = as.data.frame(t(table(srt_subset_subset$timepoint,srt_subset_subset$reannotation)))
# df = reshape2::acast(df, Var1~Var2, value.var="Freq")
# prop = apply(df,2,prop.table)
# prop = melt(prop)
# names(prop) = c("CellType","timepoint","Proportion")
# ggplot(prop, aes(fill=CellType, y=Proportion, x=timepoint)) +
#   geom_bar(position="fill", stat="identity")+
#   ggtitle("Cell Type Abundance for Non-responder")+
#   scale_fill_manual(values=c("skyblue","blue",
#                                       "purple",
#                                       "grey2","lavender",
#                                       "salmon","grey","orange","magenta","gold",
#                                       "mediumpurple1",
#                                       "limegreen","lightgreen",
#                                       "brown","pink"))

###### (2) Cell Chat ######

meta = srt@meta.data
names(meta)[4] = "samples"
meta$samples = as.factor(meta$samples)
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
index = which(meta$characteristics..response=="Responder"&meta$timepoint=="Pre")
cellchat <- createCellChat(object = as.matrix(srt@assays$RNA$data)[,index],
                           meta = meta[index,], group.by = "cell.types")
cellchat@DB <- CellChatDB
cellchat <- CellChat::subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
pdf("/Volumes/she4/hallmark/results/validation/Figure5_R_cellchat.pdf",height=10,width=10)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
# netVisual_aggregate(cellchat, signaling = "ICAM", layout = "circle")

index = which(meta$characteristics..response=="Non-responder"&meta$timepoint=="Pre")
cellchat2 <- createCellChat(object = as.matrix(srt@assays$RNA$data)[,index],
                           meta = meta[index,], group.by = "cell.types")
cellchat2@DB <- CellChatDB
cellchat2 <- CellChat::subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2, type = "triMean")
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
groupSize2 <- as.numeric(table(cellchat2@idents))
netVisual_circle(cellchat2@net$count, vertex.weight = groupSize2, weight.scale = T, label.edge= F, title.name = "Number of interactions")
pdf("/Volumes/she4/hallmark/results/validation/Figure5_NR_cellchat.pdf",height=10,width=10)
netVisual_circle(cellchat2@net$weight, vertex.weight = groupSize2, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# netVisual_aggregate(cellchat2, signaling = "ICAM", layout = "circle")
df_R = subsetCommunication(cellchat)
df_NR = subsetCommunication(cellchat2)
df_merged <- merge(df_R[,c("source","target","prob")], df_NR[,c("source","target","prob")], by = c("source", "target"), suffixes = c("_R", "_NR"))
df_merged$prob_diff <- with(df_merged, prob_R - prob_NR)
combo = expand.grid(unique(df_merged$source),unique(df_merged$target))
combo$pval = NA
for(i in 1:nrow(combo)){
  if(length(which(df_merged$source==combo[i,1]&df_merged$target==combo[i,2]))<=1){
    combo$pval[i] = 1
  }else{
    R = df_merged[which(df_merged$source==combo[i,1]&df_merged$target==combo[i,2]),3]
    NR = df_merged[which(df_merged$source==combo[i,1]&df_merged$target==combo[i,2]),4]
    combo$pval[i] = t.test(R,NR)$p.value
  }
}
combo$padj = p.adjust(combo$pval,method="fdr")
combo = combo[which(combo$padj<0.05),]
combo$R_prob = NA
combo$NR_prob = NA
for(i in 1:nrow(combo)){
  combo$R_prob[i] = sum(df_merged[which(df_merged$source==combo[i,1]&df_merged$target==combo[i,2]),3])
  combo$NR_prob[i] = sum(df_merged[which(df_merged$source==combo[i,1]&df_merged$target==combo[i,2]),4])
}
combo$diff = ifelse(combo$NR_prob-combo$R_prob>0,"NR","R")
combo$interaction = paste0(combo$Var1,"---",combo$Var2)
standardize_format <- function(input_string) {
  # Split the string by "--"
  elements <- strsplit(input_string, "---")[[1]]
  # Sort the elements alphabetically
  sorted_elements <- sort(elements,decreasing = T)
  # Join the sorted elements back together with "--"
  result <- paste(sorted_elements, collapse = "---")
  return(result)
}
for(i in 1:nrow(combo)){
  combo$interaction[i] = standardize_format(combo$interaction[i])
}
names(combo)[c(1,2,5,6,7,8)] = c("source","target","edge_R","edge_NR","enriched_in","undirected_interaction")
write.csv(combo,"/Volumes/she4/hallmark/results/validation/Figure5_differentialCellChat.csv")
pdf("/Volumes/she4/hallmark/results/validation/Figure5_differentialCellChat.pdf",height=5,width=10)
par(mfrow=c(1, 2))
g <- graph_from_data_frame(d = combo[which(combo$enriched_in=="R"),c(1,2,5)], directed = TRUE)
E(g)$width <- as.vector(scale(combo[which(combo$enriched_in=="R"),5],center=F,scale = T))
p1=plot(g,layout = layout_in_circle(g) ,edge.width = E(g)$width, layout = layout_with_fr(g), 
        vertex.size = 15,edge.arrow.size = 0.5,
        vertex.color="cornflowerblue",edge.color="black")
g <- graph_from_data_frame(d = combo[which(combo$enriched_in=="NR"),c(1,2,6)], directed = TRUE)
E(g)$width <- as.vector(scale(combo[which(combo$enriched_in=="NR"),6],center=F,scale = T))
p2=plot(g,layout = layout_in_circle(g) ,edge.width = E(g)$width, layout = layout_with_fr(g), 
        vertex.size = 15,edge.arrow.size = 0.5,
        vertex.color="cornflowerblue",edge.color="black")
dev.off()

###### (3) Detect Doublet - relate to MP ######

srt = readRDS(paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet_celltypist.rds"))
srt = subset(srt,subset = characteristics..therapy=="anti-PD1")
Idents(srt) = "cell.types"
markers = FindAllMarkers(srt,assay="RNA",only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
markers = markers[which(markers$p_val_adj<0.05),]
celltype_gs = list()
for(i in unique(markers$cluster)){
  celltype_gs[[i]] = markers[which(markers$cluster==i),"gene"]
}
ssgsea = corto::ssgsea(as.matrix(srt@assays$RNA$counts),c(gs_list,celltype_gs),minsize = 0)
ssgsea = as.data.frame(t(ssgsea))
rownames(ssgsea) = colnames(srt)
srt = AddMetaData(srt,metadata = ssgsea)
saveRDS(srt,paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet_celltypist_ALLssgsea.rds"))

srt = readRDS(paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet_celltypist_ssgsea.rds"))
srt = subset(srt,subset = characteristics..therapy=="anti-PD1")
srt_subset = subset(srt,subset=timepoint=="Pre")
doublet_srt = subset(srt_subset,subset=doublet=="doublet")
Idents(doublet_srt) = "title"

##### (3_pre) Pseudobulk #####
srt = readRDS(paste0(path_data,"/validation-data/GSE120575/seurat_umap_doublet_celltypist_ALLssgsea.rds"))
srt = subset(srt,subset = characteristics..therapy=="anti-PD1")
srt_subset = subset(srt,subset=timepoint=="Pre")
srt_sc = CreateSeuratObject(srt_subset@assays$RNA$counts,meta.data = srt_subset@meta.data)
srt_bulk = AggregateExpression(srt_sc,assays ="RNA",
                               group.by='characteristics..patinet.ID..Pre.baseline..Post..on.treatment.',
                               return.seurat = T)
ssgsea = corto::ssgsea(as.matrix(srt_bulk@assays$RNA$scale.data),gs_list,minsize = 0)
ssgsea = as.data.frame(t(ssgsea))
rownames(ssgsea) = colnames(srt_bulk)
srt_bulk = AddMetaData(srt_bulk,metadata = ssgsea)
full_data = cbind(srt_bulk@meta.data[,c("M_Eosinophil_chemotaxis",
                                        "L_CellCycle_Immune_Response",
                                        "L_Cell_Cycle_1",
                                        "L_Cell_Adhesion_2.CD52..",
                                        "M_MHC2_mediated_Lymphocyte_Activation",
                                        "L_Lipid_Localization_TCR_synapse",
                                        "L_Cytokine_Receptor_Signaling",
                                        "L_Antigen_Presentation",
                                        "L_MHC.II_mediated_immunity",
                                        "L_Cell_Adhesion_1")],
                  unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.))
full_data$response = NA
for(i in 1:nrow(full_data)){
  full_data$response[i] = unique(srt_sc$characteristics..response[which(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.==full_data$`unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)`[i])])
}
set.seed(88)
ss = umap(full_data[,1:10],n_neighbors=10)
km.final <- kmeans(full_data[,1:10], 2)
def = as.data.frame(cbind(ss$layout,full_data$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p1=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE120575 - 10 irMPs")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")

full_data = cbind(srt_bulk@meta.data[,c(61:110)],
                  unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.))
full_data$response = NA
for(i in 1:nrow(full_data)){
  full_data$response[i] = unique(srt_sc$characteristics..response[which(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.==full_data$`unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)`[i])])
}
set.seed(88)
ss = umap(full_data[,1:50],n_neighbors=10)
km.final <- kmeans(full_data[,1:50], 2)
def = as.data.frame(cbind(ss$layout,full_data$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p2=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE120575 - 50 Hallmarks")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")

full_data = cbind(srt_bulk@meta.data[,c(32:60)],
                  unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.))
full_data$response = NA
for(i in 1:nrow(full_data)){
  full_data$response[i] = unique(srt_sc$characteristics..response[which(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.==full_data$`unique(srt_sc$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)`[i])])
}
set.seed(88)
ss = umap(full_data[,1:29],n_neighbors=10)
km.final <- kmeans(full_data[,1:29], 2)
def = as.data.frame(cbind(ss$layout,full_data$response,km.final$cluster))
def$V1 = as.numeric(def$V1)
def$V2 = as.numeric(def$V2)
names(def)[3] = "Response"
p3=ggplot(def, aes(x=V1,y=V2))+
  geom_point(aes(color = Response))+
  geom_mark_ellipse(aes(color = Response), expand = unit(0.01,"mm"))+
  xlab("UMAP_1")+ylab("UMAP_2")+
  ggtitle("GSE120575 - 29 irKEGGs")+
  theme_classic()+theme(legend.position = 'bottom', legend.direction = "horizontal")

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.GSE120575_UMAP.pdf",
    height=8,width=8)
p1 + p2 + p3 + patchwork::plot_layout(guides = "collect",ncol=2) & theme(legend.position = 'bottom',legend.justification = "left")
dev.off()

f = srt_bulk@meta.data[,c("M_Eosinophil_chemotaxis",
                          "L_CellCycle_Immune_Response",
                          "L_Cell_Cycle_1",
                          "L_Cell_Adhesion_2.CD52..",
                          "M_MHC2_mediated_Lymphocyte_Activation",
                          "L_Lipid_Localization_TCR_synapse",
                          "L_Cytokine_Receptor_Signaling",
                          "L_Antigen_Presentation",
                          "L_MHC.II_mediated_immunity",
                          "L_Cell_Adhesion_1","prolif","exhaustion_full","effector","tissue_memory.full")]
testRes = corrplot::cor.mtest(f, conf.level = 0.95)
mat = cor(f)[1:10,12:14]
colnames(mat) = c("Exhaustion","Cytotoxicity","TRM")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.scRNA_baseline_pseudobulk.pdf",
    height=7,width=9)
corrplot(mat, method = 'circle', col = rev(COL2('RdBu')),
         tl.col = 'black',tl.cex = 1.5)$corrPos -> p1
text(p1$x, p1$y, round(mat, 1),cex=1.2)
dev.off()

########## (3a) Double Annotation ########
df = doublet_srt@meta.data[,c(5,9,52:62)]
df_long = melt(setDT(df),id.vars = c("title","characteristics..response"))
df_long$value = as.numeric(df_long$value)
df_long_subset = df_long[which(df_long$value>0),]
# df_long_subset = df_long
double_annot_res = matrix(0,length(unique(df_long_subset$title)),2)
for(i in 1:length(unique(df_long_subset$title))){
  double_annot = df_long_subset[which(df_long_subset$title==unique(df_long_subset$title)[i]),][order(value,decreasing = T),][c(1,2),3]
  double_annot_res[i,2] = paste0(unlist(unname(double_annot[1])),"---",unlist(unname(double_annot[2])))
}
double_annot_res[,1] = unique(df_long_subset$title)
double_annot_res = as.data.frame(double_annot_res)
double_annot_res$response = NA
for(i in 1:nrow(double_annot_res)){
  double_annot_res$response[i] = unlist(unname(df_long_subset[match(double_annot_res$V1[i],df_long_subset$title),2]))
}
standardize_format <- function(input_string) {
  
  elements <- strsplit(input_string, "---")[[1]]
  # sorted_elements <- sort(elements,decreasing =T)
  # result <- paste(sorted_elements, collapse = "---")
  # return(result)
  
  if(grepl("hausted",input_string)){
    elements <- strsplit(input_string, "---")[[1]]
    setdiff(c(1,2),grep("hausted",elements))
    result <- paste0(elements[grep("hausted",elements)], 
                    "---",
                    elements[setdiff(c(1,2),grep("hausted",elements))])
  }
  if(grepl("Memory",input_string)){
    elements <- strsplit(input_string, "---")[[1]]
    setdiff(c(1,2),grep("Memory",elements))
    result <- paste0(elements[grep("Memory",elements)], 
                     "---",
                     elements[setdiff(c(1,2),grep("Memory",elements))])
  }
  if(!grepl("hausted",input_string)&!grepl("Memory",input_string)){
    sorted_elements <- sort(elements,decreasing =T)
    result <- paste(sorted_elements, collapse = "---")
  }
  return(result)
}
for(i in 1:nrow(double_annot_res)){
  double_annot_res$V2[i] = standardize_format(double_annot_res$V2[i])
}
# doublet_of_interest = names(table(doublet_srt$doublet_annotation))[which(unname(table(doublet_srt$doublet_annotation))!=1)]
# doublet_srt = subset(doublet_srt,subset=doublet_annotation%in%doublet_of_interest)
doublet_srt$doublet_annotation = NA
for(i in 1:nrow(doublet_srt@meta.data)){
  if(length(which(double_annot_res$V1==doublet_srt$title[i]))!=0){
    doublet_srt$doublet_annotation[i] = double_annot_res[which(double_annot_res$V1==doublet_srt$title[i]),2]
  }
}
doublet_srt$doublet_annotation_group = NA
for(i in 1:nrow(doublet_srt@meta.data)){
  doublet_srt$doublet_annotation_group[i] = ifelse(grepl("Regulatory",doublet_srt$doublet_annotation[i])&grepl("Memory",doublet_srt$doublet_annotation[i]),"MemoryT-Tregs",
                                                   ifelse(grepl("Regulatory",doublet_srt$doublet_annotation[i])&!grepl("Memory",doublet_srt$doublet_annotation[i]),"Tregs-others",
                                                          ifelse(!grepl("Regulatory",doublet_srt$doublet_annotation[i])&grepl("Memory",doublet_srt$doublet_annotation[i]),"Memory-others","others-others")))
}

# doublet_srt = subset(srt_subset,subset=doublet=="doublet")
# Idents(doublet_srt) = "title"
# df = doublet_srt@meta.data[,c(5,9,52:62)]
# df_long = melt(setDT(df),id.vars = c("title","characteristics..response"))
# df_long$value = as.numeric(df_long$value)
# df_long_subset = df_long[which(df_long$value>0),]
# double_annot_res = matrix(0,length(unique(df_long_subset$title)),2)
# for(i in 1:length(unique(df_long_subset$title))){
#   double_annot = sort(df_long_subset[which(df_long_subset$title==unique(df_long_subset$title)[i]),],by=value,decreasing = T)[c(1,2),3]
#   double_annot_res[i,2] = paste0(unlist(unname(double_annot[1])),"---",unlist(unname(double_annot[2])))
# }
# double_annot_res[,1] = unique(df_long_subset$title)
# double_annot_res = as.data.frame(double_annot_res)
# double_annot_res$response = NA
# for(i in 1:nrow(double_annot_res)){
#   double_annot_res$response[i] = unlist(unname(df_long_subset[match(double_annot_res$V1[i],df_long_subset$title),2]))
# }
df = as.data.frame(t(table(double_annot_res$response,double_annot_res$V2)))
df = reshape2::acast(df, Var1~Var2, value.var="Freq")
prop = apply(df,2,prop.table)
prop = melt(prop)
names(prop) = c("Doublet","Sample","Proportion")
prop = prop[-grep("NA",prop$Doublet),]
# prop = prop[-grep("Dendritic",prop$Doublet),]
# prop = prop[-grep("Monocyte",prop$Doublet),]
# prop$Doublet = factor(prop$Doublet,
#                       levels = c("Dendritic cells---B-cells","Dendritic cells---Cytotoxicity lymphocytes",               
#                                  "Exhausted/HS CD8+ T cells---Dendritic cells","Lymphocytes exhausted/cell-cycle---B-cells",               
#                                  "Lymphocytes exhausted/cell-cycle---Exhausted CD8+ T cells","Lymphocytes---Cytotoxicity lymphocytes",                   
#                                  "Lymphocytes---Dendritic cells","Lymphocytes---Exhausted/HS CD8+ T cells",                  
#                                  "Memory T cells---B-cells","Memory T cells---Cytotoxicity lymphocytes",                
#                                  "Memory T cells---Dendritic cells","Memory T cells---Lymphocytes",    
#                                  "Regulatory T cells---Memory T cells",
#                                  "Monocytes/Macrophages---Memory T cells", 
#                                  "Monocytes/Macrophages---Cytotoxicity lymphocytes","Monocytes/Macrophages---Dendritic cells",                  
#                                  "Monocytes/Macrophages---Lymphocytes",                 
#                                  "Plasma cells---Dendritic cells","Plasma cells---Exhausted/HS CD8+ T cells",                 
#                                  "Regulatory T cells---B-cells","Regulatory T cells---Cytotoxicity lymphocytes",            
#                                  "Regulatory T cells---Dendritic cells","Regulatory T cells---Exhausted/HS CD8+ T cells",           
#                                  "Regulatory T cells---Lymphocytes","Regulatory T cells---Lymphocytes exhausted/cell-cycle",    
#                                 "Regulatory T cells---Monocytes/Macrophages" ))
pdf("/Volumes/she4/hallmark/results/validation/Figure5_doublet_enrichment.pdf",height=8,width=12)
ggplot(prop, aes(fill=Doublet, y=Proportion, x=Sample)) +
  geom_bar(position="fill", stat="identity")+
  ggtitle("Doublet annotation for each doublet at baseline")+
  scale_fill_manual(values=c("burlywood4","burlywood2",
                                         "palegreen","lightgreen","darkseagreen3","darkolivegreen1",
                                         "limegreen","darkgreen","mediumseagreen","olivedrab","seagreen3",
                                         "mistyrose","pink","plum2","plum3","hotpink","violetred3",
                                         "yellow","gold","orange",
                                         "powderblue","cornflowerblue","blue",
                                         "mediumorchid3","plum4","mediumslateblue",
                                         "mediumpurple1","blueviolet",
                                         "mediumorchid4","purple4"))
                                      # ,"burlywood4",
                                      # "antiquewhite2","aquamarine2","purple2","brown2"))
# ggplot(prop, aes(fill=Doublet, y=Proportion, x=Sample)) +
#   geom_flow(aes(alluvium = Doublet), alpha= .5)+
#   geom_col(aes(fill = Doublet), width = .35, color = "black") +
#   #geom_bar(position="fill", stat="identity", width=0.5)+
#   ggtitle("Doublet annotation for each doublet at baseline")+
#   scale_fill_manual(values=c("burlywood4","burlywood2",
#                                          "palegreen","lightgreen","darkseagreen3","darkolivegreen1","limegreen","darkgreen",
#                                          "pink","brown1","hotpink","magenta",
#                                          "yellow","gold","orange","chocolate",
#                                          "cornflowerblue","blue",
#                                          "plum2","plum3","plum4","mediumslateblue",
#                                          "mediumpurple1","blueviolet","darkmagenta",
#                                          "mediumorchid4"))+
#                                            scale_y_continuous(expand = c(0,0))
dev.off()


########## (3c) Run 3b first: Doublet enrichment R vs NR at baseline ########
df = as.data.frame(doublet_srt@meta.data[which(doublet_srt@meta.data$timepoint=="Pre"),
                                         c("orig.ident","timepoint","doublet_annotation","characteristics..response")])
df_pre = df[which(df$characteristics..response=="Responder"),]
res1 = NULL
for (i in 1:length(unique(df_pre$doublet_annotation))){
  res1 = rbind(res1,
               c(length(which(df_pre$doublet_annotation==unique(df_pre$doublet_annotation)[i])),nrow(df_pre)))
}
res1 = as.data.frame(res1)
res1$`Cell Type` = unique(df_pre$doublet_annotation)
res1$response = "Responder"

df_post = df[which(df$characteristics..response=="Non-responder"),]
res2 = NULL
for (i in 1:length(unique(df_post$doublet_annotation))){
  res2 = rbind(res2,c(length(which(df_post$doublet_annotation==unique(df_post$doublet_annotation)[i])),nrow(df_post)))
}
res2 = as.data.frame(res2)
res2$`Cell Type` = unique(df_post$doublet_annotation)
res2$response = "Non-responder"

res = as.data.frame(rbind(res1,res2))
names(res) = c("Abundance","Total","Doublet Composition","Response")
# res = res[order(res$`Doublet Composition`),]
res = res[-grep("NA",res$`Doublet Composition`),]
res = res[complete.cases(res),]
pval1 = NULL
pval2 = NULL
for (ct in 1:length(unique(res$`Doublet Composition`))){
    if(length(which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"))==0){
      dat <- data.frame(
        "celltype" = c(0,
                       res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1]),
        "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
                           sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
        row.names = c("Pre", "After"),
        stringsAsFactors = FALSE
      )
    }else if (length(which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"))==0){
      dat <- data.frame(
        "celltype" = c(res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1],
                       0),
        "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
                           sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
        row.names = c("Pre", "After"),
        stringsAsFactors = FALSE
      )
    }else{
      dat <- data.frame(
        "celltype" = c(res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1],
                       res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1]),
        "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
                           sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
        row.names = c("Pre", "After"),
        stringsAsFactors = FALSE
      )
    }
  pval1 = c(pval1,fisher.test(as.matrix(dat))$p.value)
  a = oddsratio(as.matrix(dat),log=FALSE)
  pval2 = c(pval2,1-pchisq(1/unname(coef(a)),df=1))

}
names(pval1) = unique(res$`Doublet Composition`)
pval1 = p.adjust(pval1, method = "fdr", n = length(pval1))
res$Response <- factor(res$Response, levels = c("Responder", "Non-responder"))
stat.test = as.data.frame(cbind(names(pval1),rep("Relative Abundance",13),rep("Pre",13),rep("Post",13),unname(pval1)))
names(stat.test) = c("Doublet Composition",".y.","group1","group2","p.adj")
stat.test = as_tibble(stat.test)
stat.test$`Doublet Composition`=as.factor(stat.test$`Doublet Composition`)
stat.test$p.adj.sig = ifelse(as.numeric(stat.test$p.adj)<0.00001,"****",
                             ifelse(as.numeric(stat.test$p.adj)<0.0001,"***",
                                    ifelse(as.numeric(stat.test$p.adj)<0.001,"**",
                                           ifelse(as.numeric(stat.test$p.adj)<0.05,"*","ns"))))
res$`Relative Abundance` = res$Abundance/res$Total
res$`Doublet Composition` = factor(res$`Doublet Composition`,levels = colnames(figure_df_agg)[column_order(a)])
bp <- ggbarplot(
  data = res, x = "Doublet Composition", y = "Relative Abundance",
  fill= "Response", palette = c("cornflowerblue", "pink2"),
  position = position_dodge(0.8)
)
stat.test <- stat.test %>%
  rstatix::add_xy_position(data=res,formula=`Relative Abundance`~Response, x = "dose", dodge = 0.8)
stat.test$xmin = (seq(1,52,2)+0.3)/2
stat.test$xmax = (seq(2,52,2)+0.3)/2
stat.test$p.adj = round(as.numeric(stat.test$p.adj),2)
pdf("/Volumes/she4/hallmark/results/validation/Figure5_scRNA_doublet_RvsNR_baseline.pdf",height=6,width=15)
bp +
  # stat_pvalue_manual(hide.ns = FALSE,
  #                       stat.test,  label = "p.adj.sig",
  #                       tip.length = 0.01,y.position = 0.4,bracket.size = 0.6)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
dev.off()
doublet_interested = stat.test$`Doublet Composition`
# pdf("~/desktop/doublet_annotation_MP_score.pdf",width=20,height=15,onefile = T)
# for(i in 1:length(unique(res_coef_long$mp))){
#   print(VlnPlot(doublet_srt,features = unique(res_coef_long$mp)[i],group.by = "doublet_annotation",pt.size = 0)+
#           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#           ggpubr::stat_compare_means(list=c("Memory-others","Tregs-others"),label = "p.signif")+
#           stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95))
# }
# dev.off()

########## (3b) MP activity in doublets ########
figure_df = doublet_srt@meta.data[,c(24:51,63)]
figure_df_agg = aggregate(figure_df,by=list(figure_df$doublet_annotation),FUN=mean)
figure_df_agg = figure_df_agg[-grep("NA",figure_df_agg$Group.1),]
rownames(figure_df_agg) = figure_df_agg$Group.1
figure_df_agg = figure_df_agg[,-c(1,ncol(figure_df_agg))]
# figure_df_agg = apply(figure_df_agg,2,as.numeric)
figure_df_agg = scale(figure_df_agg)
figure_df_agg = t(figure_df_agg)
figure_df_agg = figure_df_agg[unique(as.character(res_coef_long$mp)),]
figure_df_agg = figure_df_agg[,as.character(stat.test$`Doublet Composition`)]
pdf("/Volumes/she4/hallmark/results/validation/Figure5_MP_doublets.pdf",height=10,width=15)
a = Heatmap(figure_df_agg,name="MP score",show_row_dend = F,show_column_dend = F,
        cluster_rows = F,show_row_names = F,cluster_columns = T,
        column_names_max_height = unit(12, "cm"),
        row_names_max_width = unit(10, "cm"))
Heatmap(figure_df_agg,name="MP score",show_row_dend = F,show_column_dend = F,
        cluster_rows = F,show_row_names = F,cluster_columns = T,
        column_names_max_height = unit(12, "cm"),
        row_names_max_width = unit(10, "cm"))
dev.off()

# ################ GSE221553 ################
# 
# files = list.files("/Volumes/she4/hallmark/data/validation-data/GSE222446",full.names = T,pattern = ".tsv")
# mat = read.table(files[1])
# srt = CreateSeuratObject(mat)
# for(i in 2:length(files)){
#   mat = read.table(files[i])
#   srt_new = CreateSeuratObject(mat)
#   srt = merge(srt,srt_new)
#   print(i)
# }
# srt = JoinLayers(srt)
# srt$patient = NA
# for(i in 1:dim(srt)[2]){
#   srt$patient[i] = strsplit(strsplit(colnames(srt)[i],"\\.")[[1]][2],"_")[[1]][1]
# }
# srt$response = ifelse(srt$patient%in%c("patient1","patient4","patient10"),"SD",
#                       ifelse(srt$patient%in%c("patient2","patient3"),"CR",
#                              ifelse(srt$patient%in%c("patient5","patient6","patient11","patient12"),"PD","PR")))
# srt <- srt%>%
#   NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE) %>%
#   RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
# options(repr.plot.height = 2.5, repr.plot.width = 6)
# srt <- srt %>%
#   RunHarmony("patient", plot_convergence = TRUE)
# DimPlot(srt, reduction = "harmony", pt.size = .1, group.by = "patient")
# srt <- srt %>%
#   RunUMAP(reduction = "harmony", dims = 1:20) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# DimPlot(srt, reduction = "umap", group.by = "patient")
# DimPlot(srt, reduction = "umap", label = TRUE)
# FeaturePlot(srt,features = c("CD8A","CD4","FOXP3","LYZ","MS4A1","CD14","PTPRC"))
# 
# srt = subset(srt,subset=seurat_clusters%in%c(1,2,7,8,20,22,11,5)) # Get rid of malignant cells and re-cluster immune cells
# srt <- srt%>%
#   NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE) %>%
#   RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)%>%
#   FindNeighbors(dims = 1:20) %>%
#   FindClusters(resolution = 0.5) %>%
#   RunUMAP(dims = 1:20)
# DimPlot(srt, reduction = "umap", label = TRUE)
# FeaturePlot(srt,features = c("CD8A","CD4","FOXP3","LYZ","MS4A1","CD14","NKG7","GNLY"))
# saveRDS(srt,"/Volumes/she4/hallmark/data/validation-data/GSE222446/seurat_GSE221553.rds")
# writeMM(srt@assays$RNA$counts,"/Volumes/she4/hallmark/data/validation-data/GSE222446/matrix_ACT.mtx")
# 
# srt = readRDS("/rsrch4/home/bcb/she4/hallmark/data/validation-data/GSE222446/seurat_GSE221553.rds")
# celltype_gs = readRDS("/rsrch4/home/bcb/she4/hallmark/data/validation-data/GSE120575/celltype.markers.rds")
# ssgsea = corto::ssgsea(as.matrix(srt@assays$RNA$counts),c(gs_list,celltype_gs),minsize = 0)
# ssgsea = as.data.frame(t(ssgsea))
# rownames(ssgsea) = colnames(srt)
# srt = AddMetaData(srt,metadata = ssgsea)
# saveRDS(srt,paste0(path_data,"/validation-data/GSE222446/seurat_GSE221553_ALLssgsea.rds"))
# 
# srt = readRDS("/Volumes/she4/hallmark/data/validation-data/GSE222446/seurat_GSE221553_ssgsea.rds")
# doublet = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE222446/doublet_score.csv",row.names = 1)
# doublet$doublet = ifelse(doublet$X0>unname(quantile(doublet$X0, probs = 0.95)),"doublet","singlet")
# rownames(doublet) = colnames(srt)
# srt = AddMetaData(srt,metadata = doublet)
# saveRDS(srt,paste0("/Volumes/she4/hallmark/data","/validation-data/GSE222446/seurat_GSE221553_doublet_ssgsea.rds"))
# 
# doublet_srt = subset(srt,subset=doublet=="doublet")
# df = cbind(colnames(doublet_srt),doublet_srt@meta.data[,c(5,36:46)])
# names(df)[1] = "title"
# df_long = melt(setDT(df),id.vars = c("title","response"))
# df_long$value = as.numeric(df_long$value)
# df_long_subset = df_long[which(df_long$value>0),]
# # df_long_subset = df_long
# double_annot_res = matrix(0,length(unique(df_long_subset$title)),2)
# for(i in 1:length(unique(df_long_subset$title))){
#   double_annot = df_long_subset[which(df_long_subset$title==unique(df_long_subset$title)[i]),][order(value,decreasing = T),][c(1,2),3]
#   double_annot_res[i,2] = paste0(unlist(unname(double_annot[1])),"---",unlist(unname(double_annot[2])))
# }
# double_annot_res[,1] = unique(df_long_subset$title)
# double_annot_res = as.data.frame(double_annot_res)
# double_annot_res$response = NA
# for(i in 1:nrow(double_annot_res)){
#   double_annot_res$response[i] = unlist(unname(df_long_subset[match(double_annot_res$V1[i],df_long_subset$title),2]))
# }
# standardize_format <- function(input_string) {
#   # Split the string by "--"
#   elements <- strsplit(input_string, "---")[[1]]
#   # Sort the elements alphabetically
#   sorted_elements <- sort(elements,decreasing = T)
#   # Join the sorted elements back together with "--"
#   result <- paste(sorted_elements, collapse = "---")
#   return(result)
# }
# for(i in 1:nrow(double_annot_res)){
#   double_annot_res$V2[i] = standardize_format(double_annot_res$V2[i])
# }
# # doublet_of_interest = names(table(doublet_srt$doublet_annotation))[which(unname(table(doublet_srt$doublet_annotation))!=1)]
# # doublet_srt = subset(doublet_srt,subset=doublet_annotation%in%doublet_of_interest)
# double_annot_res$response_group = ifelse(double_annot_res$response%in%c("PD","SD"),"Non-responder","Responder")
# doublet_srt$response_group = ifelse(doublet_srt$response%in%c("PD","SD"),"Non-responder","Responder")
# doublet_srt$doublet_annotation = NA
# for(i in 1:nrow(doublet_srt@meta.data)){
#   if(length(which(double_annot_res$V1==colnames(doublet_srt)[i]))!=0){
#     doublet_srt$doublet_annotation[i] = double_annot_res[which(double_annot_res$V1==colnames(doublet_srt)[i]),2]
#   }
# }
# doublet_srt$doublet_annotation_group = NA
# for(i in 1:nrow(doublet_srt@meta.data)){
#   doublet_srt$doublet_annotation_group[i] = ifelse(grepl("Regulatory",doublet_srt$doublet_annotation[i])&grepl("Memory",doublet_srt$doublet_annotation[i]),"MemoryT-Tregs",
#                                                    ifelse(grepl("Regulatory",doublet_srt$doublet_annotation[i])&!grepl("Memory",doublet_srt$doublet_annotation[i]),"Tregs-others",
#                                                           ifelse(!grepl("Regulatory",doublet_srt$doublet_annotation[i])&grepl("Memory",doublet_srt$doublet_annotation[i]),"Memory-others","others-others")))
# }
# 
# doublet_interested
# 
# df = as.data.frame(doublet_srt@meta.data[,c("orig.ident","doublet_annotation","response_group")])
# df_pre = df[which(df$response_group=="Responder"),]
# res1 = NULL
# for (i in 1:length(unique(df_pre$doublet_annotation))){
#   res1 = rbind(res1,
#                c(length(which(df_pre$doublet_annotation==unique(df_pre$doublet_annotation)[i])),nrow(df_pre)))
# }
# res1 = as.data.frame(res1)
# res1$`Cell Type` = unique(df_pre$doublet_annotation)
# res1$response = "Responder"
# 
# df_post = df[which(df$response_group=="Non-responder"),]
# res2 = NULL
# for (i in 1:length(unique(df_post$doublet_annotation))){
#   res2 = rbind(res2,c(length(which(df_post$doublet_annotation==unique(df_post$doublet_annotation)[i])),nrow(df_post)))
# }
# res2 = as.data.frame(res2)
# res2$`Cell Type` = unique(df_post$doublet_annotation)
# res2$response = "Non-responder"
# 
# res = as.data.frame(rbind(res1,res2))
# names(res) = c("Abundance","Total","Doublet Composition","Response")
# res = res[which(res$`Doublet Composition` %in% doublet_interested),]
# # res = res[order(res$`Doublet Composition`),]
# # res = res[-grep("NA",res$`Doublet Composition`),]
# res = res[complete.cases(res),]
# pval1 = NULL
# pval2 = NULL
# for (ct in 1:length(unique(res$`Doublet Composition`))){
#   if(length(which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"))==0){
#     dat <- data.frame(
#       "celltype" = c(0,
#                      res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1]),
#       "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
#                          sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
#       row.names = c("Pre", "After"),
#       stringsAsFactors = FALSE
#     )
#   }else if (length(which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"))==0){
#     dat <- data.frame(
#       "celltype" = c(res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1],
#                      0),
#       "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
#                          sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
#       row.names = c("Pre", "After"),
#       stringsAsFactors = FALSE
#     )
#   }else{
#     dat <- data.frame(
#       "celltype" = c(res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1],
#                      res[which(res$`Doublet Composition`==unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1]),
#       "non_celltype" = c(sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Responder"),1]),
#                          sum(res[which(res$`Doublet Composition`!=unique(res$`Doublet Composition`)[ct]&res$Response=="Non-responder"),1])),
#       row.names = c("Pre", "After"),
#       stringsAsFactors = FALSE
#     )
#   }
#   pval1 = c(pval1,fisher.test(as.matrix(dat))$p.value)
#   a = oddsratio(as.matrix(dat),log=FALSE)
#   pval2 = c(pval2,1-pchisq(1/unname(coef(a)),df=1))
# 
# }
# names(pval1) = unique(res$`Doublet Composition`)
# pval1 = p.adjust(pval1, method = "fdr", n = length(pval1))
# res$Response <- factor(res$Response, levels = c( "Responder","Non-responder"))
# 
# stat.test = as.data.frame(cbind(names(pval1),rep("Relative Abundance",25),rep("Pre",25),rep("Post",25),unname(pval1)))
# names(stat.test) = c("Doublet Composition",".y.","group1","group2","p.adj")
# stat.test = as_tibble(stat.test)
# stat.test$`Doublet Composition`=as.factor(stat.test$`Doublet Composition`)
# stat.test$p.adj.sig = ifelse(as.numeric(stat.test$p.adj)<0.00001,"****",
#                              ifelse(as.numeric(stat.test$p.adj)<0.0001,"***",
#                                     ifelse(as.numeric(stat.test$p.adj)<0.001,"**",
#                                            ifelse(as.numeric(stat.test$p.adj)<0.05,"*","ns"))))
# res$`Relative Abundance` = res$Abundance/res$Total
# res$`Doublet Composition` = factor(res$`Doublet Composition`,levels = colnames(figure_df_agg)[column_order(a)])
# res = as.data.frame(rbind(as.matrix(res),rbind(c(0,941,"Plasma cells---Exhausted/HS CD8+ T cells","Responder",0),
#                       c(0,459,"Plasma cells---Exhausted/HS CD8+ T cells","Non-responder",0))))
# res$`Relative Abundance` = as.numeric(res$`Relative Abundance`)
# res$`Doublet Composition` = factor(res$`Doublet Composition`,levels = colnames(figure_df_agg)[column_order(a)])
# # stat.test = stat.test[match(colnames(figure_df_agg)[column_order(a)],stat.test$`Doublet Composition`),]
# # stat.test = stat.test[-c(1,25),]
# bp <- ggbarplot(
#   data = res, x = "Doublet Composition", y = "Relative Abundance",
#   fill= "Response", palette = c("cornflowerblue","pink2"),
#   position = position_dodge(0.8)
# )
# stat.test <- stat.test %>%
#   rstatix::add_xy_position(data=res,formula=`Relative Abundance`~Response, x = "dose", dodge = 0.8)
# stat.test$xmin = (seq(1,48,2)+0.3)/2
# stat.test$xmax = (seq(2,48,2)+0.3)/2
# stat.test$p.adj = round(as.numeric(stat.test$p.adj),2)
# pdf("/Volumes/she4/hallmark/results/validation/Figure5_ACTscRNA_doublet_RvsNR_baseline.pdf",height=6,width=15)
# bp +
#   # stat_pvalue_manual(hide.ns = FALSE,
#   #                       stat.test,  label = "p.adj.sig",
#   #                       tip.length = 0.01,y.position = 0.4,bracket.size = 0.6)+
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
#   theme_minimal()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
# dev.off()
