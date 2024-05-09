path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

library(dplyr)
library(AnnotationDbi)
library(ComplexHeatmap)
library(org.Hs.eg.db)

# Our MPs
annotation_L = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res.csv"))
annotation_M = read.csv(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_res.csv"))

gs_list_M = list()
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names)){
  gs_list_M[[i]] = read.csv(file_names[i],row.names = 1)[,1]
}
names(gs_list_M) = names(annotation_M)

gs_list_L = list()
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = glob2rx("scaled*_based.csv")) #where you have your files
for (i in 1:length(file_names)){
  gs_list_L[[i]] = read.csv(file_names[i],row.names = 1)[,1]
}
names(gs_list_L) = names(annotation_L)

# Spectra gene sets
spectra = read.csv("/Volumes/she4/hallmark/data/spectra_genesets.csv")
spectra = spectra[,c(1,6)]
spectra_geneset = list()
for(i in 1:nrow(spectra)){
  a = gsub("[[]","",spectra$gene_set[i])
  a = gsub("[]]","",a)
  a = gsub("[']","",a)
  spectra_geneset[[i]] = unlist(strsplit(a,", "))
}
names(spectra_geneset) = spectra$gene_set_name

# Hallmark
gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
H_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(H_list) = unique(gene_sets$gs_name)

# C7
gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "C7"))
C7_list = list()
for (i in 1:length(unique(gene_sets$gs_name))){
  C7_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
}
names(C7_list) = unique(gene_sets$gs_name)

# KEGG
kegg_sets = read.csv(paste0(path_data,"kegg_immune_sets.csv"))[,-1]
KEGG_list = list()
for (i in 1:length(unique(kegg_sets$name))){
  KEGG_list[[i]] = unique(kegg_sets[which(kegg_sets$name==unique(kegg_sets$name)[i]),"symbol"])
}
names(KEGG_list) = unique(kegg_sets$name)

# WIKI
# wp_immune_pathway_ID = c("WP5038","WP5078","WP4329","WP4705","WP5039","WP4559","WP5116","WP4582","WP4542","WP3863","WP2328","WP4197","WP4449","WP4754",
#                          "WP1836","WP2794","WP2719","WP4460","WP3577","WP1799","WP2671","WP1835","WP2694","WP4998","WP4484","WP1798","WP2688","WP2763",
#                          "WP2551","WP585","WP4584","WP3865","WP5130","WP5225","WP3945","WP3937","WP1895","WP3858","WP3890","WP4914","WP195","WP5055","WP4146","WP4091",
#                          "WP5092","WP5198","WP5049","WP5010","WP3591","WP3287","WP4853","WP3869","WP4891","WP4298","WP4747","WP4217","WP4630","WP3594","WP4847","WP5222",
#                          "WP4666","WP3802","WP15","WP5115")
# 
# gmt <- 'https://wikipathways-data.wmcloud.org/20220510/gmt/wikipathways-20220510-gmt-Homo_sapiens.gmt'
# read.gmt.wp <- function(gmtfile) {
#   clusterProfiler::read.gmt(gmtfile) %>%
#     tidyr::separate(.data$term, c("name","version","wpid","species"), "%")
# }
# wp <- read.gmt.wp(gmt)
# wp = wp[which(wp$wpid %in% wp_immune_pathway_ID),] # only 41/64 were found
# cols <- c("SYMBOL", "GENENAME","ENTREZID")
# map = AnnotationDbi::select(org.Hs.eg.db, 
#                             keys=wp$gene, 
#                             columns=cols, keytype="ENTREZID")
# wp$symbol = map[which(map$ENTREZID%in%wp$gene),"SYMBOL"]
# WP_list = list()
# for (i in 1:length(unique(wp$wpid))){
#   WP_list[[i]] = unique(wp[which(wp$wpid==unique(wp$wpid)[i]),"symbol"])
# }
# names(WP_list) = unique(wp$name)
# saveRDS(WP_list,paste0("/Volumes/she4/hallmark/data/wp.rds"))
WP_list = readRDS("/rsrch4/home/bcb/she4/hallmark/data/wp.rds")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
jaccard7 <- function(a, b) {
  intersection = length(intersect(a, b))
  return (intersection)
}

spectra = read.csv("/Volumes/she4/hallmark/data/spectra_genesets.csv")
spectra$celltype_major = factor(spectra$celltype_major,
                                levels=c("all-cells","DC","Endothelial","ILC","Macrophage","Mast","NK","B","T"))
spectra = spectra[order(spectra$celltype_major),]
table(spectra$celltype_major)
spectra_geneset = spectra_geneset[spectra$gene_set_name]

## ------------------------------------------ between MP and spectra ------------------------------------------ ##
res = matrix(0,9,length(spectra_geneset))
pvalue = matrix(0,9,length(spectra_geneset))
for (i in 1:length(gs_list_M)){
  for (j in 1:231){
    res[i,j] = jaccard(gs_list_M[[i]], spectra_geneset[[j]])
    iter = 0
    ja = NULL
    while (iter < 100){
      iter = iter + 1
      ja = c(ja,jaccard(gs_list_M[[i]], sample(unique(unlist(spectra_geneset)),length(spectra_geneset[[j]]))))
      print(iter)
    }
    pvalue[i,j] = ifelse(jaccard(gs_list_M[[i]], spectra_geneset[[j]])<0.05,1,
                         sum(ja>jaccard(gs_list_M[[i]], spectra_geneset[[j]]))/100)
  }
}
pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=9,ncol=231)
colnames(res) = colnames(pvalue) = names(spectra_geneset)
rownames(res) = rownames(pvalue) = names(gs_list_M)
# pdf(paste0(path_result,"new_overlap_spectra_M_density.pdf"),height = 2.4,width = 5)
# den = as.data.frame(as.vector(res))
# names(den)="value"
# ggplot(den,aes(x=value))+geom_density(color="darkblue", fill="lightblue")+
#   xlab("Jaccard distance between spectra gene sets and myeloid MPs")+ylab("Frequency")
# dev.off()
ann_col <- data.frame(spectra$celltype_major)
colnames(ann_col) <- c("Cell Type")
col = c("lightblue2","lightcoral","gold","steelblue","pink","green3","lavender","purple","orange")
names(col) = unique(spectra$celltype_major)
col = list(col)
names(col) = c("Cell Type")
colAnn <- HeatmapAnnotation(annotation_width = unit(c(1, 4), 'cm'),
                            foo = anno_block(gp = gpar(fill = c("lightblue2","lightcoral","gold","steelblue","pink","green3","lavender","purple","orange")), 
                                             labels = unique(spectra$celltype_major)),
                            gap = unit(1, 'mm'),show_annotation_name=F)
pdf(paste0(path_result,"new_overlap_spectra_M.pdf"),height = 8,width = 30)
Heatmap(res,name="Jaccard Distance",
        column_title = "Overlap between Myeloid MPs and spectra pathways",
        show_row_names = T,show_column_names = T,
        show_row_dend = F,show_column_dend = F,
        cluster_columns = F,cluster_rows = F,
        column_names_max_height = unit(9,"cm"),
        row_names_max_width = unit(10.5,"cm"),
        top_annotation = colAnn,
        column_split=factor(spectra$celltype_major),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

res = matrix(0,19,length(names(spectra_geneset)))
pvalue = matrix(0,19,length(names(spectra_geneset)))
for (i in 1:length(gs_list_L)){
  for (j in 1:231){
    res[i,j] = jaccard(gs_list_L[[i]], spectra_geneset[[j]])
    iter = 0
    ja = NULL
    while (iter < 100){
      iter = iter + 1
      ja = c(ja,jaccard(gs_list_L[[i]], sample(unique(unlist(spectra_geneset)),length(spectra_geneset[[j]]))))
      print(iter)
    }
    pvalue[i,j] = ifelse(jaccard(gs_list_L[[i]], spectra_geneset[[j]])<0.05,1,
                         sum(ja>jaccard(gs_list_L[[i]], spectra_geneset[[j]]))/100)
  }
}
pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),
                nrow=19,ncol=231)
colnames(res) = colnames(pvalue) = names(spectra_geneset)
rownames(res) = rownames(pvalue) = names(gs_list_L)
# pdf(paste0(path_result,"new_overlap_spectra_L_density.pdf"),height = 2.4,width = 5)
# den = as.data.frame(as.vector(res))
# names(den)="value"
# ggplot(den,aes(x=value))+geom_density(color="darkblue", fill="lightblue")+
#   xlab("Jaccard distance between spectra gene sets and lymphoid MPs")+ylab("Frequency")
# dev.off()
ann_col <- data.frame(spectra$celltype_major)
colnames(ann_col) <- c("Cell Type")
col = c("lightblue2","lightcoral","gold","steelblue","pink","green3","lavender","purple","orange")
names(col) = unique(spectra$celltype_major)
col = list(col)
names(col) = c("Cell Type")
colAnn <- HeatmapAnnotation(annotation_width = unit(c(1, 4), 'cm'),
                            foo = anno_block(gp = gpar(fill = c("lightblue2","lightcoral","gold","steelblue","pink","green3","lavender","purple","orange")), 
                                             labels = unique(spectra$celltype_major)),
                            gap = unit(1, 'mm'),show_annotation_name=F)

pdf(paste0(path_result,"new_overlap_spectra_L.pdf"),height = 12,width = 35)
Heatmap(res,name="Jaccard Distance",
        column_title = "Overlap between Lymphoid MPs and spectra pathways",
        show_row_names = T,show_column_names = T,
        show_row_dend = F,show_column_dend = F,
        cluster_columns = F,cluster_rows = F,
        column_names_max_height = unit(9,"cm"),
        row_names_max_width = unit(10.5,"cm"),
        top_annotation = colAnn,
        column_split=factor(spectra$celltype_major),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

## ------------------------------------------ between C7 and C7 ------------------------------------------ ##
res = matrix(0,5219,5219)
for (i in 1:length(C7_list)){
  for (j in 1:5219){
    res[i,j] = jaccard(C7_list[[i]], C7_list[[j]])
  }
}
colnames(res) = names(C7_list)
rownames(res) = names(C7_list)
write.csv(res,"~/desktop/between_C7.csv")

# col = c("white","skyblue","blue")
# pdf(paste0(path_result,"overlap_itself_C7.pdf"),height = 8,width = 10)
# diag(res[1:2000,1:2000]) = 0
# Heatmap(res[1:2000,1:2000],name="Jaccard Distance",col = col,
#         column_title = "Overlap across Lymphoid MP pathways",
#         show_row_names =F,show_column_names = F,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = T,
#         column_names_max_height = unit(9,"cm"))
# dev.off()

## ------------------------------------------ between MP and C7 ------------------------------------------ ##
# pvalue: probabilty of observing something as extreme as the the observed one: 
# = the probability of having higher overlap. If < 0.05, then the overlap is significant, should be reported
# res = matrix(0,9,length(C7_list))
# pvalue = matrix(0,9,length(C7_list))

# res = matrix(0,9,length(C7_list))
# pvalue = matrix(0,9,length(C7_list))
# for (i in 1:length(gs_list_M)){
#   for (j in x:y){ # 1:length(C7_list)
#     res[i,j] = jaccard(gs_list_M[[i]], C7_list[[j]])
#     iter = 0
#     ja = NULL
#     while (iter < 100){
#       iter = iter + 1
#       ja = c(ja,jaccard(gs_list_M[[i]], sample(unique(unlist(C7_list)),length(C7_list[[j]]))))
#       print(iter)
#     }
#     pvalue[i,j] = ifelse(jaccard(gs_list_M[[i]], C7_list[[j]])<0.02,1,
#                          formatC(sum(ja>jaccard(gs_list_M[[i]], C7_list[[j]]))/100, format = "e"))
#   }
# }
# #pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=9,ncol=29)
# colnames(res) = names(C7_list)
# rownames(res) = names(gs_list_M)
# write.csv(res,paste0("/rsrch4/home/bcb/she4/hallmark/results_overlap/M_",x,"_",y,"_res.csv"))
# write.csv(pvalue,paste0("/rsrch4/home/bcb/she4/hallmark/results_overlap/M_",x,"_",y,"_pvalue.csv"))

# pdf(paste0(path_result,"new_overlap_C7_M.pdf"),height = 8,width = 23)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Myeloid MPs and ImmuneSigDB pathways",
#         show_row_names = T,show_column_names = F,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(9,"cm"),
#         row_names_max_width = unit(10.5,"cm"))
# dev.off()

res = matrix(0,19,length(C7_list))
pvalue = matrix(0,19,length(C7_list))
for (i in 1:length(gs_list_L)){
  for (j in x:y){ # 1:length(C7_list)
    res[i,j] = jaccard(gs_list_L[[i]], C7_list[[j]])
    iter = 0
    ja = NULL
    while (iter < 100){
      iter = iter + 1
      ja = c(ja,jaccard(gs_list_L[[i]], sample(unique(unlist(C7_list)),length(C7_list[[j]]))))
    }
    pvalue[i,j] = ifelse(jaccard(gs_list_L[[i]], C7_list[[j]])<0.02,1,
                         formatC(sum(ja>jaccard(gs_list_L[[i]], C7_list[[j]]))/100, format = "e"))
  }
}
colnames(res) = names(C7_list)
rownames(res) = names(gs_list_L)
write.csv(res,paste0("/rsrch4/home/bcb/she4/hallmark/results_overlap/L_",x,"_",y,"_res.csv"))
write.csv(pvalue,paste0("/rsrch4/home/bcb/she4/hallmark/results_overlap/L_",x,"_",y,"_pvalue.csv"))

#pdf(paste0(path_result,"new_overlap_C7_L.pdf"),height = 8,width = 23)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Lymphoid MPs and ImmuneSigDB pathways",
#         show_row_names = T,show_column_names = F,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(9,"cm"),
#         row_names_max_width = unit(10.5,"cm"))
# dev.off()

# ## ------------------------------------------ between MP and MP ------------------------------------------ ##
# 
# res = matrix(0,19,19)
# for (i in 1:length(gs_list_L)){
#   for (j in 1:19){
#     res[i,j] = jaccard(gs_list_L[[i]], gs_list_L[[j]])
#   }
# }
# 
# # res = read.csv(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/MPs.csv"))[,-1]
# # overlap_matrix = outer(1:ncol(res), 1:ncol(res), FUN = Vectorize(function(a, b)(
# #   length(intersect(unlist(res[,a]), unlist(res[,b]))))))
# colnames(res) = names(gs_list_L)
# rownames(res) = names(gs_list_L)
# col = c("white","skyblue","blue")
# pdf(paste0(path_result,"overlap_itself_L.pdf"),height = 8,width = 10)
# Heatmap(res,name="Jaccard Distance",col = col,
#         column_title = "Overlap across Lymphoid MP pathways",
#         show_row_names =F,show_column_names = F,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = T,
#         column_names_max_height = unit(9,"cm"),
#         heatmap_legend_param = list(title = "Overlap", at = c(0,0.5,1), column_title_gp = gpar(fontsize = 15),
#                                     row_title_gp = gpar(fontsize = 10), col = col),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#                                       grid.text(sprintf("%.2f", res[i, j]), x, y, gp = gpar(fontsize = 15))
#                                     })
# dev.off()
# 
# res = matrix(0,9,9)
# for (i in 1:length(gs_list_M)){
#   for (j in 1:9){
#     res[i,j] = jaccard(gs_list_M[[i]], gs_list_M[[j]])
#   }
# }
# colnames(res) = names(gs_list_M)
# rownames(res) = names(gs_list_M)
# col = c("white","skyblue","blue")
# pdf(paste0(path_result,"overlap_itself_M.pdf"),height = 8,width = 10)
# Heatmap(res,name="Jaccard Distance",col = col,
#         column_title = "Overlap across Myeloid MP pathways",
#         show_row_names = F,show_column_names = F,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = T,
#         column_names_max_height = unit(9,"cm"),
#         heatmap_legend_param = list(title = "Overlap", at = c(0,0.5,1), column_title_gp = gpar(fontsize = 15),
#                                     row_title_gp = gpar(fontsize = 10), col = col),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", res[i, j]), x, y, gp = gpar(fontsize = 15))
#         })
# dev.off()
# 
# 
## ------------------------------------------ between MP and hallmark ------------------------------------------ ##
res = matrix(0,9,length(H_list))
pvalue = matrix(0,9,length(H_list))
for (i in 1:length(gs_list_M)){
  for (j in 1:50){
    res[i,j] = jaccard(gs_list_M[[i]], H_list[[j]])
    iter = 0
    ja = NULL
    while (iter < 100){
      iter = iter + 1
      ja = c(ja,jaccard(gs_list_M[[i]], sample(unique(unlist(H_list)),length(H_list[[j]]))))
      print(iter)
    }
    pvalue[i,j] = ifelse(jaccard(gs_list_M[[i]], H_list[[j]])<0.05,1,
                         sum(ja>jaccard(gs_list_M[[i]], H_list[[j]]))/100)
  }
}
pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=9,ncol=50)
colnames(res) = colnames(pvalue) = substr(names(H_list),10,nchar(names(H_list)))
rownames(res) = rownames(pvalue) = names(gs_list_M)
pdf(paste0(path_result,"new_overlap_H_M.pdf"),height = 8,width = 14)
Heatmap(res,name="Jaccard Distance",
        column_title = "Overlap between Myeloid MPs and Hallmark pathways",
        show_row_names = T,show_column_names = T,
        show_row_dend = F,show_column_dend = F,
        cluster_columns = F,cluster_rows = F,
        column_names_max_height = unit(9,"cm"),
        row_names_max_width = unit(10.5,"cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

res = matrix(0,19,length(H_list))
pvalue = matrix(0,19,length(H_list))
for (i in 1:length(gs_list_L)){
  for (j in 1:50){
    res[i,j] = jaccard(gs_list_L[[i]], H_list[[j]])
    iter = 0
    ja = NULL
    while (iter < 100){
      iter = iter + 1
      ja = c(ja,jaccard(gs_list_L[[i]], sample(unique(unlist(H_list)),length(H_list[[j]]))))
    }
    pvalue[i,j] = ifelse(jaccard(gs_list_L[[i]], H_list[[j]])<0.05,1,
                         sum(ja>jaccard(gs_list_L[[i]], H_list[[j]]))/100)
  }
}
pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=19,ncol=50)
colnames(res) = colnames(pvalue) = substr(names(H_list),10,nchar(names(H_list)))
rownames(res) = rownames(pvalue) = names(gs_list_L)
pdf(paste0(path_result,"new_overlap_H_L.pdf"),height = 8,width = 14)
Heatmap(res,name="Jaccard Distance",
        column_title = "Overlap between Lymphoid MPs and Hallmark pathways",
        show_row_names = T,show_column_names = T,
        show_row_dend = F,show_column_dend = F,
        cluster_columns = T,cluster_rows = F,
        column_names_max_height = unit(9,"cm"),
        row_names_max_width = unit(10.5,"cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
        })
dev.off()

# ## ------------------------------------------ between MP and kegg   ------------------------------------------ ##
# res = matrix(0,9,length(KEGG_list))
# pvalue = matrix(0,9,length(KEGG_list))
# for (i in 1:length(gs_list_M)){
#   for (j in 1:29){
#     res[i,j] = jaccard(gs_list_M[[i]], KEGG_list[[j]])
#     iter = 0
#     ja = NULL
#     while (iter < 100){
#       iter = iter + 1
#       ja = c(ja,jaccard(gs_list_M[[i]], sample(unique(unlist(KEGG_list)),length(KEGG_list[[j]]))))
#     }
#     pvalue[i,j] = ifelse(jaccard(gs_list_M[[i]], KEGG_list[[j]])<0.02,1,
#                          sum(ja>jaccard(gs_list_M[[i]], KEGG_list[[j]]))/100)
#   }
# }
# pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=9,ncol=29)
# colnames(res) = colnames(pvalue) = names(KEGG_list)
# rownames(res) = rownames(pvalue) = names(gs_list_M)
# pdf(paste0(path_result,"new_overlap_KEGG_M.pdf"),height = 8,width = 14)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Myeloid MPs and KEGG Immune pathways",
#         show_row_names = T,show_column_names = T,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(9,"cm"),
#         row_names_max_width = unit(10.5,"cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# 
# res = matrix(0,19,length(KEGG_list))
# pvalue = matrix(0,19,length(KEGG_list))
# for (i in 1:length(gs_list_L)){
#   for (j in 1:29){
#     res[i,j] = jaccard(gs_list_L[[i]], KEGG_list[[j]])
#     iter = 0
#     ja = NULL
#     while (iter < 100){
#       iter = iter + 1
#       ja = c(ja,jaccard(gs_list_L[[i]], sample(unique(unlist(KEGG_list)),length(KEGG_list[[j]]))))
#     }
#     pvalue[i,j] = ifelse(jaccard(gs_list_L[[i]], KEGG_list[[j]])<0.02,1,
#                          sum(ja>jaccard(gs_list_L[[i]], KEGG_list[[j]]))/100)
#   }
# }
# pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=19,ncol=29)
# colnames(res) = colnames(pvalue) = names(KEGG_list)
# rownames(res) = rownames(pvalue) = names(gs_list_L)
# pdf(paste0(path_result,"new_overlap_KEGG_L.pdf"),height = 8,width = 14)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Lymphoid MPs and KEGG Immune pathways",
#         show_row_names = T,show_column_names = T,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(9,"cm"),
#         row_names_max_width = unit(10.5,"cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# 
# ## ------------------------------------------ between MP and wikipath ------------------------------------------ ##
# res = matrix(0,9,length(WP_list))
# pvalue = matrix(0,9,length(WP_list))
# for (i in 1:length(gs_list_M)){
#   for (j in 1:41){
#     res[i,j] = jaccard(gs_list_M[[i]], WP_list[[j]])
#     iter = 0
#     ja = NULL
#     while (iter < 100){
#       iter = iter + 1
#       ja = c(ja,jaccard(gs_list_M[[i]], sample(unique(unlist(WP_list)),length(WP_list[[j]]))))
#     }
#     pvalue[i,j] = ifelse(jaccard(gs_list_M[[i]], WP_list[[j]])<0.02,1,
#                          sum(ja>jaccard(gs_list_M[[i]], WP_list[[j]]))/100)
#   }
# }
# pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=9,ncol=41)
# colnames(res) = colnames(pvalue) = names(WP_list)
# rownames(res) = rownames(pvalue) = names(gs_list_M)
# pdf(paste0(path_result,"new_overlap_wiki_M.pdf"),height = 11,width = 14)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Myeloid MPs and Wiki Immune pathways",
#         show_row_names = T,show_column_names = T,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(15,"cm"),
#         row_names_max_width = unit(10.5,"cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()
# 
# res = matrix(0,19,length(WP_list))
# pvalue = matrix(0,19,length(WP_list))
# for (i in 1:length(gs_list_L)){
#   for (j in 1:41){
#     res[i,j] = jaccard(gs_list_L[[i]], WP_list[[j]])
#     iter = 0
#     ja = NULL
#     while (iter < 100){
#       iter = iter + 1
#       ja = c(ja,jaccard(gs_list_L[[i]], sample(unique(unlist(WP_list)),length(WP_list[[j]]))))
#     }
#     pvalue[i,j] = ifelse(jaccard(gs_list_L[[i]], WP_list[[j]])<0.02,1,
#                          sum(ja>jaccard(gs_list_L[[i]], WP_list[[j]]))/100)
#   }
# }
# pvalue = matrix(p.adjust(as.vector(pvalue),method="fdr",n=length(as.vector(pvalue))),nrow=19,ncol=41)
# colnames(res) = colnames(pvalue) = names(WP_list)
# rownames(res) = rownames(pvalue) = names(gs_list_L)
# pdf(paste0(path_result,"new_overlap_wiki_L.pdf"),height = 11,width = 14)
# Heatmap(res,name="Jaccard Distance",
#         column_title = "Overlap between Lymphoid MPs and Wiki Immune pathways",
#         show_row_names = T,show_column_names = T,
#         show_row_dend = F,show_column_dend = F,
#         cluster_columns = T,cluster_rows = F,
#         column_names_max_height = unit(15,"cm"),
#         row_names_max_width = unit(10.5,"cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(pvalue[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         })
# dev.off()




