library(data.table)
library(ComplexHeatmap)

a = list.files("/Volumes/she4/hallmark/results_overlap",pattern = "pvalue",full.names = T)
a = a[grep("L_",a)]
b = list.files("/Volumes/she4/hallmark/results_overlap",pattern = "res",full.names = T)
b = b[grep("L_",b)]

pvalue_all = matrix(NA,19,5219)
for (i in 1:length(a)){
  start = as.numeric(strsplit(a[i],"_")[[1]][3])
  end = as.numeric(strsplit(a[i],"_")[[1]][4])-1
  pvalue_all[,c(start:end)] = as.matrix(read.csv(a[i],row.names = 1)[,c(start:end)])
  print(i)
}

# names(C7_list)[which(apply(pvalue_all,2,sum)!=9)]

res_all = matrix(NA,19,5219)
for (i in 1:length(b)){
  start = as.numeric(strsplit(b[i],"_")[[1]][3])
  end = as.numeric(strsplit(b[i],"_")[[1]][4])
  fillin = as.matrix(fread(b[i],select = c(start:end)))
  fillin = as.data.frame(fillin)
  fillin = fillin[,-1]
  fillin = as.matrix(fillin)
  fillin = apply(fillin,2,as.numeric)
  #fillin = read.csv(b[i],row.names = 1,col.names = paste0("V",c(start:end)))
  res_all[,c(start:(end-1))] = fillin
  print(i)
}
View(res_all[,which(apply(pvalue_all,2,sum)!=19)])
colnames(res_all) = names(C7_list)
rownames(res_all) = names(gs_list_L)

# Two requirements (1) pvalue<0.05 (2) overlap>0.08
sig_mat = res_all[,which(apply(pvalue_all,2,sum)!=19)]
colnames(sig_mat) = colnames(res_all)[which(apply(pvalue_all,2,sum)!=19)]
sig_big_mat = sig_mat[,unique(which(sig_mat>0.06, arr.ind = TRUE)[,2])]

sig_p = pvalue_all[,which(apply(pvalue_all,2,sum)!=19)]
colnames(sig_p) = colnames(res_all)[which(apply(pvalue_all,2,sum)!=19)]
sig_big_p = sig_p[,unique(which(sig_mat>0.06, arr.ind = TRUE)[,2])]

pdf(paste0("~/desktop/","new_overlap_C7_L.pdf"),height = 18,width = 21)
Heatmap(t(sig_big_mat),name="Jaccard Distance",
        column_title = "Overlap between Lymphoid MPs and ImmuneSigDB pathways",
        show_row_names = T,show_column_names = T,
        show_row_dend = F,show_column_dend = F,
        cluster_columns = T,cluster_rows = T,
        column_names_max_height = unit(10.5,"cm"),
        row_names_max_width = unit(35,"cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(t(sig_big_p)[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 15))
          })
dev.off()


