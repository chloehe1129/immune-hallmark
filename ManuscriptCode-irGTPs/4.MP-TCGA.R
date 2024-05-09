#path_code = "/Volumes/she4/hallmark/code/"
path_code = "/rsrch4/home/bcb/she4/hallmark/code/"

#path_data = "/Volumes/she4/hallmark/data/"
path_data = "/rsrch4/home/bcb/she4/hallmark/data/"

#path_result = "/Volumes/she4/hallmark/results/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

library(readxl)
library(AnnotationDbi)
library(corto)
library(org.Hs.eg.db)
library(survival)
# library(GSVA)
library(survminer)
library(ComplexHeatmap)
library(glmnet)
library(survival)

# -------------------------------- Input the TCGA matrix -------------------------------- #

# # Expression matrix
# count = read.table(paste0(path_data,"validation-data/TCGA/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"))
# names(count) = count[1,]
# names(count) = substr(names(count),1,12)
# count = count[-1,]
# 
# a = unlist(strsplit(count[,1], '|', fixed = TRUE))
# count[,1] = a[seq(2,length(a),by=2)]
# count = count[!duplicated(count[,1]),]
# rownames(count) = count[,1]
# count = count[,-1]
# 
# # Metadata
# metadata = read_excel(paste0(path_data,"/validation-data/TCGA/TCGA-CDR-SupplementalTableS1.xlsx"),
#                  sheet="TCGA-CDR")[,-1]
# names(metadata)[3]="age"
# metadata = metadata[which(metadata$bcr_patient_barcode%in%names(count)),]
# write.csv(metadata,paste0(path_data,"/validation-data/TCGA/metadata.csv"))
# 
# # Expression matrix
# cols <- c("SYMBOL", "GENENAME","ENTREZID")
# map = AnnotationDbi::select(org.Hs.eg.db, 
#                             keys=rownames(count), 
#                             columns=cols, keytype="ENTREZID")
# map = map[!duplicated(map$ENTREZID),]
# count$symbol = map$SYMBOL
# count = count[!duplicated(count$symbol),]
# count = count[!is.na(count$symbol),]
# rownames(count) = count$symbol
# count = count[,-which(names(count)=="symbol")]
# count = count[,which(names(count)%in%metadata$bcr_patient_barcode)]
# write.csv(count,paste0(path_data,"/validation-data/TCGA/pancancer_count.csv"))

# -------------------------------- Create MP score -------------------------------- #

metadata = read.csv(paste0(path_data,"/validation-data/TCGA/metadata.csv"))[,-1]
count = read.csv(paste0(path_data,"/validation-data/TCGA/pancancer_count.csv"))
rownames(count) = count[,1]
count = count[,-1]

# Our MPs
gs_list = list()
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = "_based") #where you have your files
for (i in 1:length(file_names)){
  gs_list[[i]] = read.csv(file_names[i])[,-1]
}
names(gs_list) = paste0("MP",seq(1,13,1))
# 
# # Hallmark
# gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
# H_list = list()
# for (i in 1:length(unique(gene_sets$gs_name))){
#   H_list[[i]] = unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[i]),"gene_symbol"])
# }
# names(H_list) = unique(gene_sets$gs_name)
# 
# # Kegg 
# kegg_sets = read.csv(paste0("/rsrch4/home/bcb/she4/hallmark/data/","kegg_immune_sets.csv"))[,-1]
# KEGG_list = list()
# for (i in 1:length(unique(kegg_sets$name))){
#   KEGG_list[[i]] = unique(kegg_sets[which(kegg_sets$name==unique(kegg_sets$name)[i]),"symbol"])
# }
# names(KEGG_list) = unique(kegg_sets$name)
# 
# ssGSEA - mp score
mat = apply(count,2,as.numeric)
rownames(mat) = rownames(count)
ssGSEA_MP = corto::ssgsea(mat,gs_list)
write.csv(as.data.frame(t(ssGSEA_MP)),paste0(path_result,"humanMyeloid_TCGA_ssgsea_MP.csv"))

# 
# # ssGSEA - kegg score
# ssGSEA_KEGG = corto::ssgsea(mat, KEGG_list) 
# write.csv(as.data.frame(t(ssGSEA_KEGG)),paste0(path_result,"TCGA_ssgsea_KEGG.csv"))
# 
# # ssGSEA - hallmark score
# ssGSEA_H = corto::ssgsea(mat, H_list)
# write.csv(as.data.frame(t(ssGSEA_H)),paste0(path_result,"TCGA_ssgsea_hallmark.csv"))

# -------------------------------- Run survival LASSP analysis -------------------------------- #

ssGSEA_MP = read.csv(paste0(path_result,"humanMyeloid_TCGA_ssgsea_MP.csv"))[,-1]
df = cbind(ssGSEA_MP,metadata)
result_MP = matrix(0,length(unique(df$type)),13*2)

res = NULL
for (cancer in 1:length(unique(df$type))){
  p = NULL    
  subset = df[which(df$type==unique(df$type)[cancer]),]
  subset = subset[,c(1:13,16,19,38,39)]
  subset = subset[complete.cases(subset),]

  if(sum(subset$OS.time<=0)!=0){subset = subset[-which(subset$OS.time==0),]}
  survobj = Surv(subset$OS.time,subset$OS)
  x = subset[,c(1:(ncol(subset)-2))] # MPs,age,stage
  x$ajcc_pathologic_tumor_stage = as.numeric(factor(x$ajcc_pathologic_tumor_stage))
  
  cv.fit <- cv.glmnet(as.matrix(x), y=survobj, family = "cox")
  fit <- glmnet(x, y=survobj, family = "cox", lambda=cv.fit$lambda.min)
  num = length(rownames(coef(fit))[-which(coef(fit)==0)])
  res = rbind(res,
              cbind(rep(unique(df$type)[cancer],num),rownames(coef(fit))[-which(coef(fit)==0)],round(as.numeric(coef(fit)[-which(coef(fit)==0)]),3)))
  
  # pdf(paste0(path_result,unique(df$type)[cancer],"_MP_kmplot.pdf"),height = 25,width=20,onefile=T)
  # print(arrange_ggsurvplots(p[1:20],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[13:40],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[41:60],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[61:80],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[81:100],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[101:120],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[113:127],ncol = 4,nrow = 5))
  # dev.off()
}
write.csv(res, paste0(path_result,"humanMyeloid_TCGA_MP_survival_LASSO.csv"))
myeloid_lasso = read.csv("~/desktop/myeloid.csv",row.names = 1)
lymphoid_lasso = read.csv("~/desktop/lymphoid.csv",row.names = 1)

# pdf("~/desktop/myeloid_lasso.pdf")
# Heatmap(as.matrix(myeloid_lasso),na_col = "gray",cluster_rows = F,cluster_columns = F,
#         name = "LASSO Coef for Myeloid MPs")
# dev.off()
# pdf("~/desktop/lymphoid_lasso.pdf")
# Heatmap(as.matrix(lymphoid_lasso),na_col = "gray",cluster_rows = F,cluster_columns = F,
#         name = "LASSO Coef for Lymphoid MPs")
# dev.off()

# -------------------------------- Run survival analysis -------------------------------- #

ssGSEA_MP = read.csv(paste0(path_result,"humanMyeloid_TCGA_ssgsea_MP.csv"))[,-1]
df = cbind(ssGSEA_MP,metadata)
result_MP = matrix(0,length(unique(df$type)),21*2)
for (cancer in 1:length(unique(df$type))){
  p = NULL    
  subset = df[which(df$type==unique(df$type)[cancer]),]
  survobj = Surv(subset$OS.time,subset$OS)
  
  for (i in 1:21){
    
    #subset$MP_group = as.factor(ifelse(subset[,i]>unname(summary(subset[,i])[3]),"High","Low"))
    #subset$MP_group = relevel(subset$MP_group, ref = "Low")
    # new_df <- with(subset,data.frame(MP_group = c("Low","High")))
    # new_df$MP_group =as.factor(new_df$MP_group) 
    # new_df$MP_group = relevel(new_df$MP_group, ref = "Low")
    if(length(unique(subset$ajcc_pathologic_tumor_stage))==1){
      cox_fit_MP <- coxph(survobj~subset[,i]+age, data=subset)
    }else{
      cox_fit_MP <- coxph(survobj~subset[,i]+age+ajcc_pathologic_tumor_stage, data=subset)
    }
    #fit_MP <- survfit(cox_fit_MP, newdata = new_df)
    #p[[i]] = ggsurvplot(fit_MP,data=new_df, conf.int = TRUE, legend.labs=c(paste0("Low MP",i," expression"), paste0("High MP",i," expression")),
    #           ggtheme = theme_minimal(),title = paste0("TCGA Cancer: ",unique(df$type)[cancer]),
    #           pval=paste0(round(summary(cox_fit_MP)$coef[1,5],3),ifelse(round(summary(cox_fit_MP)$coef[1,5],3)<0.05,"**","")))
    result_MP[cancer,c((i*2)-1,i*2)] = c(summary(cox_fit_MP)$coef[1,2],summary(cox_fit_MP)$coef[1,5])
  }
  # pdf(paste0(path_result,unique(df$type)[cancer],"_MP_kmplot.pdf"),height = 25,width=20,onefile=T)
  # print(arrange_ggsurvplots(p[1:20],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[13:40],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[41:60],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[61:80],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[81:100],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[101:120],ncol = 4,nrow = 5))
  # print(arrange_ggsurvplots(p[113:127],ncol = 4,nrow = 5))
  # dev.off()
}
rownames(result_MP) = unique(df$type)
write.csv(result_MP, paste0(path_result,"humanMyeloid_TCGA_MP_survival.csv"))


# ssGSEA_H = read.csv(paste0(path_result,"TCGA_ssgsea_hallmark.csv"))[,-1]
# df = cbind(ssGSEA_H,metadata)
# result_H = matrix(0,length(unique(df$type)),50*2)
# for (cancer in 1:length(unique(df$type))){
#   p = NULL
#   subset = df[which(df$type==unique(df$type)[cancer]),]
#   survobj = Surv(subset$OS.time,subset$OS)
#   
#   for (i in 1:50){
#     subset$H_group = as.factor(ifelse(subset[,i]>unname(summary(subset[,i])[3]),"High","Low"))
#     subset$H_group = relevel(subset$H_group, ref = "Low")
#     # new_df <- with(subset,data.frame(H_group = c("Low","High")))
#     # new_df$H_group =as.factor(new_df$H_group) 
#     # new_df$H_group = relevel(new_df$H_group, ref = "Low")
#     cox_fit_H <- coxph(survobj~H_group+age, data=subset)
#     #fit_H <- survfit(cox_fit_H, newdata = new_df)
#     #p[[i]] = ggsurvplot(fit_H,data=new_df, conf.int = TRUE, legend.labs=c(paste0("Low ",names(subset)[i]," expression"), paste0("High ",names(subset)[i]," expression")),
#     #                         ggtheme = theme_minimal(),title = paste0("TCGA Cancer: ",unique(df$type)[cancer]),
#     #                         pval=paste0(round(summary(cox_fit_H)$coef[1,5],3),ifelse(round(summary(cox_fit_H)$coef[1,5],3)<0.05,"**","")))
#     result_H[cancer,c((i*2)-1,i*2)] = c(summary(cox_fit_H)$coef[1,2],summary(cox_fit_H)$coef[1,5])
#   }
#   # pdf(paste0(path_result,unique(df$type)[cancer],"_H_kmplot.pdf"),height = 25,width=20,onefile=T)
#   # print(arrange_ggsurvplots(p[1:20],ncol = 4,nrow = 5))
#   # print(arrange_ggsurvplots(p[13:40],ncol = 4,nrow = 5))
#   # print(arrange_ggsurvplots(p[41:50],ncol = 4,nrow = 5))
#   # dev.off()
# }
# rownames(result_H) = unique(df$type)
# write.csv(result_H, paste0(path_result,"TCGA_H_survival1.csv"))
# 
# 
# ssGSEA_KEGG = read.csv(paste0(path_result,"TCGA_ssgsea_KEGG.csv"))[,-1]
# df = cbind(ssGSEA_KEGG,metadata)
# result_KEGG = matrix(0,length(unique(df$type)),29*2)
# for (cancer in 1:length(unique(df$type))){
#   p = NULL
#   subset = df[which(df$type==unique(df$type)[cancer]),]
#   survobj = Surv(subset$OS.time,subset$OS)
#   
#   for (i in 1:29){
#     subset$KEGG_group = as.factor(ifelse(subset[,i]>unname(summary(subset[,i])[3]),"High","Low"))
#     subset$KEGG_group = relevel(subset$KEGG_group, ref = "Low")
#     # new_df <- with(subset,data.frame(KEGG_group = c("Low","High")))
#     # new_df$KEGG_group =as.factor(new_df$KEGG_group) 
#     # new_df$KEGG_group = relevel(new_df$KEGG_group, ref = "Low")
#     cox_fit_KEGG <- coxph(survobj~KEGG_group+age, data=subset)
#     #fit_KEGG <- survfit(cox_fit_KEGG, newdata = new_df)
#     #p[[i]] = ggsurvplot(fit_KEGG,data=new_df, conf.int = TRUE, legend.labs=c(paste0("Low ",names(subset)[i]," expression"), paste0("High ",names(subset)[i]," expression")),
#     #                         ggtheme = theme_minimal(),title = paste0("TCGA Cancer: ",unique(df$type)[cancer]),
#     #                         pval=paste0(round(summary(cox_fit_KEGG)$coef[1,5],3),ifelse(round(summary(cox_fit_KEGG)$coef[1,5],3)<0.05,"**","")))
#     result_KEGG[cancer,c((i*2)-1,i*2)] = c(summary(cox_fit_KEGG)$coef[1,2],summary(cox_fit_KEGG)$coef[1,5])
#   }
#   # pdf(paste0(path_result,unique(df$type)[cancer],"_KEGG_kmplot.pdf"),height = 25,width=25,onefile=F)
#   # print(arrange_ggsurvplots(p,ncol = 5,nrow = 6))
#   # dev.off()
# }
# rownames(result_KEGG) = unique(df$type)
# write.csv(result_KEGG, paste0(path_result,"TCGA_KEGG_survival1.csv"))

# -------------------------------- Plot survival analysis results -------------------------------- #
# HR = result_KEGG[,c(seq(1,58,2))]
# rownames(HR)=rownames(result_KEGG)
# colnames(HR) = names(ssGSEA_KEGG)
# pval = result_KEGG[,c(seq(2,58,2))]
# rownames(pval)=rownames(result_KEGG)
# colnames(pval) = names(ssGSEA_KEGG)
# pdf(paste0(path_result,"TCGA_KEGG_survival1.pdf"),height =12,width=16,onefile=T)
# Heatmap(HR,name="Hazard Ratio",
#       cluster_rows = F,cluster_columns = T,
#       show_row_names = T,show_column_names = T,show_row_dend = F,show_column_dend = T,
#       column_names_max_height = unit(12, "cm"),
#       cell_fun = function(j, i, x, y, width, height, fill) {
#         grid.text(ifelse(pval[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 15))
#       },
#       column_title = "KEGG pathway association with TCGA cancer-specific survival")
# dev.off()

result_MP = read.csv(paste0(path_result,"humanMyeloid_TCGA_MP_survival.csv"),row.names = 1)
HR = result_MP[,c(seq(1,21*2,2))]
HR_1 = HR
rownames(HR_1) = rownames(result_MP)
colnames(HR_1) = paste0("MP",seq(1,21,1))
pval = result_MP[,c(seq(2,21*2,2))]
pval_1 = t(apply(pval,1,p.adjust,method="fdr"))
#pval_1 = matrix(p.adjust(as.vector(pval_1),method="fdr"),nrow=nrow(pval_1),ncol = ncol(pval_1))
rownames(pval_1)=rownames(result_MP)
colnames(pval_1) = paste0("MP",seq(1,13,1))
pdf(paste0(path_result,"humanMyeloid_TCGA_MP_survival1.pdf"),height =12,width=16,onefile=T)
Heatmap(as.matrix(HR_1),name="Hazard Ratio",
        cluster_rows = F,cluster_columns = F,
        show_row_names = T,show_column_names = T,show_row_dend = F,show_column_dend = T,
        column_names_max_height = unit(12, "cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(pval_1[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 15))
        },
        column_title = "MP pathway association with TCGA cancer-specific survival")
dev.off()

# HR = result_H[,c(seq(1,50*2,2))]
# rownames(HR)=rownames(result_H)
# colnames(HR) = names(ssGSEA_H)
# pval = result_H[,c(seq(2,50*2,2))]
# rownames(pval)=rownames(result_H)
# colnames(pval) = names(ssGSEA_H)
# pdf(paste0(path_result,"TCGA_H_survival1.pdf"),height =12,width=16,onefile=T)
# Heatmap(HR,name="Hazard Ratio",
#         cluster_rows = F,cluster_columns = T,
#         show_row_names = T,show_column_names = T,show_row_dend = F,show_column_dend = T,
#         column_names_max_height = unit(12, "cm"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i, j]<0.05,"*",""), x, y, gp = gpar(fontsize = 15))
#         },
#         column_title = "Hallmark pathway association with TCGA cancer-specific survival")
# dev.off()

# -------------------------------- Association between MPs and Immune abundance -------------------------------- #

metadata = read.csv(paste0(path_data,"/validation-data/TCGA/metadata.csv"))[,-1]
metadata$TCGA.Participant.Barcode = metadata$bcr_patient_barcode
ssGSEA_MP = read.csv(paste0(path_result,"humanMyeloid_TCGA_ssgsea_MP.csv"))[,-1]
rownames(ssGSEA_MP) = metadata$bcr_patient_barcode
ssGSEA_MP$TCGA.Participant.Barcode = rownames(ssGSEA_MP)
cluster = read.csv(paste0(path_data,"validation-data/TCGA/1-s2.0-S1074761318301213-mmc2.csv"))
cluster$TCGA.Participant.Barcode = cluster$TCGA.Participant.Barcode
# ssGSEA_MP = ssGSEA_MP[which(rownames(ssGSEA_MP)%in%cluster$TCGA.Participant.Barcode),]
# ssGSEA_MP$TCGA.Participant.Barcode = rownames(ssGSEA_MP)

### Correlation between cell type abundance
# df_immune = df[,c("Macrophages.M2","T.Cells.CD4.Memory.Activated","NK.Cells.Activated",
#                   "B.Cells.Naive","Th2.Cells","T.Cells.gamma.delta","T.Cells.CD4.Memory.Resting",
#                   "Th1.Cells","Macrophages.M1","Th17.Cells","Dendritic.Cells.Activated",
#                   "T.Cells.Follicular.Helper","T.Cells.CD4.Naive","B.Cells.Memory","Mast.Cells.Activated",
#                   "Neutrophils","NK.Cells.Resting","T.Cells.Regulatory.Tregs","Plasma.Cells","Eosinophils",
#                   "Mast.Cells.Resting","Monocytes","T.Cells.CD8","Dendritic.Cells.Resting")]
# df_immune = df_immune[complete.cases(df_immune),]
###

df = merge(ssGSEA_MP,merge(cluster,metadata, by="TCGA.Participant.Barcode"),by="TCGA.Participant.Barcode")
df$hot_cold = as.factor(ifelse(df$Immune.Subtype=="C4"|df$Immune.Subtype=="C5","0","1")) # cold = 0, hot = 1
df = df[!is.na(df$hot_cold),]

# 1. MP and Cell type abundance association
variables = c("Th1.Cells","Th2.Cells","Th17.Cells",
              "T.Cells.CD4.Memory.Activated","T.Cells.CD4.Memory.Resting",
              "T.Cells.CD4.Naive","T.Cells.CD8",
              "T.Cells.Follicular.Helper","T.Cells.gamma.delta","T.Cells.Regulatory.Tregs",
              "NK.Cells.Activated","NK.Cells.Resting",
              "B.Cells.Memory","B.Cells.Naive","Plasma.Cells",
              "Dendritic.Cells.Activated","Dendritic.Cells.Resting",
              "Macrophages.M1","Macrophages.M2","Mast.Cells.Activated","Mast.Cells.Resting",
              "Monocytes","Eosinophils","Neutrophils")
for (v in 1:length(variables)){
  mat = matrix(0,13,4)
  for (i in 1:13){
    res = lm(df[,variables[v]]~df[,i+1]+age+type,data=df)
    mat[i,] = unname(coef(summary(res))[2,])
  }
  rownames(mat) = paste0("MP",seq(1,13,1))
  mat = as.data.frame(mat)
  names(mat) = c("Estimate","Std.Error","Z.value","p.value")
  mat$fdr.pval = p.adjust(mat$p.value, method = "fdr", n = length(mat$p.value))
  
  assign(paste0("result_",variables[v]),mat)
  write.csv(mat,paste0(path_result,"humanMyeloid_TCGA_MP_association_",variables[v],".csv"))
}
mat = as.matrix(cbind(result_Th1.Cells[,1],result_Th17.Cells[,1],result_Th2.Cells[,1],
                      result_T.Cells.CD4.Memory.Activated[,1],result_T.Cells.CD4.Memory.Resting[,1],
                      result_T.Cells.CD4.Naive[,1],result_T.Cells.CD8[,1],
                      result_T.Cells.Follicular.Helper[,1],result_T.Cells.gamma.delta[,1],result_T.Cells.Regulatory.Tregs[,1],
                      result_NK.Cells.Activated[,1],result_NK.Cells.Resting[,1],
                      result_B.Cells.Memory[,1],result_B.Cells.Naive[,1],result_Plasma.Cells[,1],
                      result_Dendritic.Cells.Activated[,1],result_Dendritic.Cells.Resting[,1],
                      result_Macrophages.M1[,1],result_Macrophages.M2[,1],
                      result_Mast.Cells.Activated[,1],result_Mast.Cells.Resting[,1],
                      result_Monocytes[,1],result_Eosinophils[,1],result_Neutrophils[,1]))
pval = as.matrix(cbind(result_Th1.Cells[,4],result_Th17.Cells[,4],result_Th2.Cells[,4],
                       result_T.Cells.CD4.Memory.Activated[,4],result_T.Cells.CD4.Memory.Resting[,4],
                       result_T.Cells.CD4.Naive[,4],result_T.Cells.CD8[,4],
                       result_T.Cells.Follicular.Helper[,4],result_T.Cells.gamma.delta[,4],result_T.Cells.Regulatory.Tregs[,4],
                       result_NK.Cells.Activated[,4],result_NK.Cells.Resting[,4],
                       result_B.Cells.Memory[,4],result_B.Cells.Naive[,4],result_Plasma.Cells[,4],
                       result_Dendritic.Cells.Activated[,4],result_Dendritic.Cells.Resting[,4],
                       result_Macrophages.M1[,4],result_Macrophages.M2[,4],
                       result_Mast.Cells.Activated[,4],result_Mast.Cells.Resting[,4],
                       result_Monocytes[,4],result_Eosinophils[,4],result_Neutrophils[,4]))
pval = matrix(p.adjust(as.vector(pval),method="fdr"),nrow=nrow(pval),ncol = ncol(pval))
rownames(mat) = rownames(pval) = paste0("MP",seq(1,13,1))
colnames(mat) = colnames(pval) = c("Th1.Cells","Th2.Cells","Th17.Cells",
                  "T.Cells.CD4.Memory.Activated","T.Cells.CD4.Memory.Resting",
                  "T.Cells.CD4.Naive","T.Cells.CD8",
                  "T.Cells.Follicular.Helper","T.Cells.gamma.delta","T.Cells.Regulatory.Tregs",
                  "NK.Cells.Activated","NK.Cells.Resting",
                  "B.Cells.Memory","B.Cells.Naive","Plasma.Cells",
                  "Dendritic.Cells.Activated","Dendritic.Cells.Resting",
                  "Macrophages.M1","Macrophages.M2","Mast.Cells.Activated","Mast.Cells.Resting",
                  "Monocytes","Eosinophils","Neutrophils")
mat2 = apply(mat,2,norm<-function(x){return ((x - mean(x)) / sd(x))})

pdf(paste0(path_result,"humanMyeloid_TCGA_MP_association_w_celltype.pdf"),width=16,height=19)
set.seed(1129)
Heatmap(mat2,name="Coefficient",
        cluster_rows = F,
        cluster_columns = T,
        column_title = "Association between MP and cell type abundance",
        show_row_names = T,show_row_dend = T,
        show_column_names = T,
        column_names_rot=90,
        cell_fun = function(j,i, x, y, width, height, fill) {
          grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 20))
        },
        column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
        column_title_gp = gpar(fontsize = 20,fontface="bold"),
        row_names_max_width = unit(16, "cm"),
        column_names_max_height = unit(19, "cm"))
dev.off()

# 2. MPs and others (continous outcome)
variables = c("Lymphocyte.Infiltration.Signature.Score","Intratumor.Heterogeneity","IFN.gamma.Response",
              "TGF.beta.Response","Aneuploidy.Score","Nonsilent.Mutation.Rate")
for (v in 1:length(variables)){
  mat = matrix(0,13,4)
  for (i in 1:13){
    res = lm(df[,variables[v]]~df[,i+1]+age+type,data=df)
    mat[i,] = unname(coef(summary(res))[2,])
  }
  rownames(mat) = paste0("MP",seq(13,1))
  mat = as.data.frame(mat)
  names(mat) = c("Estimate","Std.Error","Z.value","p.value")
  mat$fdr.pval = p.adjust(mat$p.value, method = "fdr", n = length(mat$p.value))
  
  assign(paste0("result_",variables[v]),mat)
  write.csv(mat,paste0(path_result,"humanMyeloid_TCGA_MP_association_",variables[v],".csv"))
}

# 3. MPs and Hot/Cold
result_HC = matrix(0,13,4)
df$hot_cold <- relevel(df$hot_cold, ref = "0")
for (i in 1:13){
  res = glm(hot_cold~df[,i+1]+age+type,data=df,family = "binomial")
  result_HC[i,] = unname(coef(summary(res))[2,])
}
rownames(result_HC) = paste0("MP",seq(1,13,1))
result_HC = as.data.frame(result_HC)
names(result_HC) = c("Estimate","Std.Error","Z.value","p.value")
result_HC$fdr.pval = p.adjust(result_HC$p.value, method = "fdr", n = length(result_HC$p.value))
write.csv(result_HC,paste0(path_result,"humanMyeloid_TCGA_MP_association_ColdHot.csv"))

mat = as.matrix(cbind(result_Lymphocyte.Infiltration.Signature.Score[,1],result_Intratumor.Heterogeneity[,1],
                      result_IFN.gamma.Response[,1],result_TGF.beta.Response[,1],result_Aneuploidy.Score[,1],
                      result_Nonsilent.Mutation.Rate[,1],result_HC[,1]))
pval = as.matrix(cbind(result_Lymphocyte.Infiltration.Signature.Score[,4],result_Intratumor.Heterogeneity[,4],
                       result_IFN.gamma.Response[,4],result_TGF.beta.Response[,4],result_Aneuploidy.Score[,4],
                       result_Nonsilent.Mutation.Rate[,4],result_HC[,4]))
pval = matrix(p.adjust(as.vector(pval),method="fdr"),nrow=nrow(pval),ncol = ncol(pval))
rownames(mat) = rownames(pval) = paste0("MP",seq(1,13,1))
colnames(mat) = colnames(pval) = c("Lymphocyte.Infiltration.Signature.Score","Intratumor.Heterogeneity","IFN.gamma.Response" ,
                                   "TGF.beta.Response","Aneuploidy.Score","Nonsilent.Mutation.Rate",
                                   "Hot_Cold")
mat2 = apply(mat,2,norm<-function(x){return ((x - mean(x)) / sd(x))})

pdf(paste0(path_result,"humanMyeloid_TCGA_MP_association_w_immuneScores.pdf"),width=8,height=15)
set.seed(1129)
Heatmap(mat2,name="Coefficient",
        cluster_rows = F,
        cluster_columns = T,
        column_title = "Association between MP and Immune pathways",
        show_row_names = T,show_row_dend = T,
        show_column_names = T,
        column_names_rot=90,
        cell_fun = function(j,i, x, y, width, height, fill) {
          grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 20))
        },
        column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
        column_title_gp = gpar(fontsize = 20,fontface="bold"),
        row_names_max_width = unit(16, "cm"),
        column_names_max_height = unit(19, "cm"))
dev.off()


############################### BELOW ARE NOT NEEDED ANYMORE ##############################

# # 4. MPs and OS
# result_OS = matrix(0,127,4)
# survobj = Surv(df$OS.time,df$OS.y)
# for (i in 1:127){
#   res = coxph(survobj~df[,i+1]+age+type,data=df)
#   result_OS[i,] = unname(coef(summary(res))[2,c(1,3:5)])
# }
# rownames(result_OS) = paste0("MP",seq(1,127,1))
# result_OS = as.data.frame(result_OS)
# names(result_OS) = c("Estimate","Std.Error","Z.value","p.value")
# result_OS$fdr.pval = p.adjust(result_OS$p.value, method = "fdr", n = length(result_OS$p.value))
# write.csv(result_OS,paste0(path_result,"humanMyeloid_TCGA_MP_association_OS.csv"))
# 
# # 5. MPs and SNV.Neoantigens
# result_SNV = matrix(0,127,4)
# for (i in 1:127){
#   res = glm(SNV.Neoantigens~df[,i+1]+age+type,data=df,family = "poisson")
#   result_SNV[i,] = unname(coef(summary(res))[2,])
# }
# rownames(result_SNV) = paste0("MP",seq(1,127,1))
# result_SNV = as.data.frame(result_SNV)
# names(result_SNV) = c("Estimate","Std.Error","Z.value","p.value")
# result_SNV$fdr.pval = p.adjust(result_SNV$p.value, method = "fdr", n = length(result_SNV$p.value))
# write.csv(result_SNV,paste0(path_result,"humanMyeloid_TCGA_MP_association_SNV.csv"))
# 
# # 6. MPs and Indel.Neoantigens
# result_indel = matrix(0,127,4)
# for (i in 1:127){
#   res = glm(Indel.Neoantigens~df[,i+1]+age+type,data=df,family = "poisson")
#   result_indel[i,] = unname(coef(summary(res))[2,])
# }
# rownames(result_indel) = paste0("MP",seq(1,127,1))
# result_indel = as.data.frame(result_indel)
# names(result_indel) = c("Estimate","Std.Error","Z.value","p.value")
# result_indel$fdr.pval = p.adjust(result_indel$p.value, method = "fdr", n = length(result_indel$p.value))
# write.csv(result_indel,paste0(path_result,"humanMyeloid_TCGA_MP_association_indel.csv"))
# 
# 
# mat = as.matrix(cbind(result_Lymphocyte.Infiltration.Signature.Score[,1],result_Intratumor.Heterogeneity[,1],
#                       result_IFN.gamma.Response[,1],result_TGF.beta.Response[,1],result_Aneuploidy.Score[,1],
#                       result_Nonsilent.Mutation.Rate[,1],
#                       result_HC[,4],result_OS[,1],result_SNV[,1],result_indel[,1]))
# pval = as.matrix(cbind(result_Lymphocyte.Infiltration.Signature.Score[,4],result_Intratumor.Heterogeneity[,4],
#                        result_IFN.gamma.Response[,4],result_TGF.beta.Response[,4],result_Aneuploidy.Score[,4],
#                        result_Nonsilent.Mutation.Rate[,4],
#                        result_HC[,4],result_OS[,1],result_SNV[,4],result_indel[,4]))
# rownames(mat) = rownames(pval) = paste0("MP",seq(1,127,1))
# colnames(mat) = colnames(pval) = c("Lymphocyte.Infiltration.Signature.Score","Intratumor.Heterogeneity","IFN.gamma.Response" ,
#                                    "TGF.beta.Response","Aneuploidy.Score","Nonsilent.Mutation.Rate",
#                                    "Hot_Cold","Hazard Rate","SNV.Neoantigens","Indel.Neoantigens")
# mat2 = apply(mat,2,norm<-function(x){return ((x - mean(x)) / sd(x))})
#     
# # factoextra::fviz_nbclust(mat2, kmeans, method = "wss")
# pdf(paste0(path_result,"humanMyeloid_TCGA_MP_association_w_pathway.pdf"),width=16,height=19)
# set.seed(1129)
# Heatmap(mat2,name="Coefficient",
#         cluster_rows = T,
#         cluster_columns = T,
#         column_title = "Association between MP and Immune pathways",
#         show_row_names = T,show_row_dend = T,
#         show_column_names = T,
#         column_names_rot=90,
#         cell_fun = function(j,i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         },
#         column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
#         column_title_gp = gpar(fontsize = 20,fontface="bold"),
#         row_names_max_width = unit(16, "cm"),
#         column_names_max_height = unit(19, "cm"))
# dev.off()
# 
# mat = as.matrix(cbind(result_Th1.Cells[,1],result_Th17.Cells[,1],result_Th2.Cells[,1],
#                       result_T.Cells.CD4.Memory.Activated[,1],result_T.Cells.CD4.Memory.Resting[,1],
#                       result_T.Cells.CD4.Naive[,1],result_T.Cells.CD8[,1],
#                       result_T.Cells.Follicular.Helper[,1],result_T.Cells.gamma.delta[,1],result_T.Cells.Regulatory.Tregs[,1],
#                       result_NK.Cells.Activated[,1],result_NK.Cells.Resting[,1],
#                       result_B.Cells.Memory[,1],result_B.Cells.Naive[,1],
#                       result_Lymphocyte.Infiltration.Signature.Score[,1],
#                       result_Intratumor.Heterogeneity[,1],
#                       result_HC[,1],result_OS[,1]))
# pval = as.matrix(cbind(result_Th1.Cells[,4],result_Th17.Cells[,4],result_Th2.Cells[,4],
#                        result_T.Cells.CD4.Memory.Activated[,4],result_T.Cells.CD4.Memory.Resting[,4],
#                        result_T.Cells.CD4.Naive[,4],result_T.Cells.CD8[,4],
#                        result_T.Cells.Follicular.Helper[,4],result_T.Cells.gamma.delta[,4],result_T.Cells.Regulatory.Tregs[,4],
#                        result_NK.Cells.Activated[,4],result_NK.Cells.Resting[,4],
#                        result_B.Cells.Memory[,4],result_B.Cells.Naive[,4],
#                        result_Lymphocyte.Infiltration.Signature.Score[,4],
#                        result_Intratumor.Heterogeneity[,4],
#                        result_HC[,4],result_OS[,4]))
# rownames(mat) = rownames(pval) = paste0("MP",seq(1,127,1))
# colnames(mat) = colnames(pval) = c("Th1.Cells","Th2.Cells","Th17.Cells",
#                                    "T.Cells.CD4.Memory.Activated","T.Cells.CD4.Memory.Resting",
#                                    "T.Cells.CD4.Naive","T.Cells.CD8",
#                                    "T.Cells.Follicular.Helper","T.Cells.gamma.delta","T.Cells.Regulatory.Tregs",
#                                    "NK.Cells.Activated","NK.Cells.Resting","B.Cells.Memory",
#                                    "B.Cells.Naive","Lymphocyte.Infiltration.Signature.Score",
#                                    "Intratumor.Heterogeneity","Hot_Cold","Hazard Rate")
# mat2 = apply(mat,2,norm<-function(x){return ((x - mean(x)) / sd(x))})
# 
# # factoextra::fviz_nbclust(mat2, kmeans, method = "wss")
# pdf(paste0(path_result,"humanMyeloid_TCGA_MP_association_w_lymphocyte.pdf"),width=16,height=19)
# set.seed(1129)
# Heatmap(mat2,name="Coefficient",
#         cluster_rows = T,
#         cluster_columns = T,
#         column_title = "Association between MP and lymphocyte",
#         show_row_names = T,show_row_dend = T,
#         show_column_names = T,
#         column_names_rot=90,
#         cell_fun = function(j,i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         },
#         column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
#         column_title_gp = gpar(fontsize = 20,fontface="bold"),
#         row_names_max_width = unit(16, "cm"),
#         column_names_max_height = unit(19, "cm"))
# dev.off()
# 
# 
# # -------------------------------- Study MP Clusters -------------------------------- #
# #write.csv(mat2,"~/desktop/mat2.csv")
# set.seed(1129)
# h=Heatmap(mat2,name="Coefficient",
#         cluster_rows = T,
#         cluster_columns = T,
#         column_title = "Association between MP and Immune Response",
#         show_row_names = T,show_row_dend = T,
#         show_column_names = T,
#         column_names_rot=30,
#         cell_fun = function(j,i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         },
#         row_km=8,row_gap = unit(4, "mm"),
#         column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
#         column_title_gp = gpar(fontsize = 20,fontface="bold"),
#         row_names_max_width = unit(16, "cm"),
#         cluster_row_slices = FALSE)
# 
# set.seed(1129)
# cluster = row_order(h)
# data = as.data.frame(cbind(paste0("MP",unlist(cluster)),c(rep("MP_C1",length(cluster[[1]])),
#                                      rep("MP_C2",length(cluster[[2]])),
#                                      rep("MP_C3",length(cluster[[3]])),
#                                      rep("MP_C4",length(cluster[[4]])),
#                                      rep("MP_C5",length(cluster[[5]])),
#                                      rep("MP_C6",length(cluster[[6]])),
#                                      rep("MP_C7",length(cluster[[7]])),
#                                      rep("MP_C8",length(cluster[[8]])))))
# write.csv(data,paste0(path_result,"MP_clusters.csv"))
# 
# MC = list()
# for (i in 1:8){
#   MPC1 = data[which(data$V2==unique(data$V2)[i]),"V1"]
#   res = NULL
#   for (l in 1:length(MPC1)){
#     res = rbind(res,read_excel(paste0(path_result,"c7_NMF_programs/cluster_methodTirosh2/res.xlsx"),sheet=MPC1[l])[,c(1:4)])
#   }
#   MC[[i]] = head(unlist(unname(res[order(res$Counts,decreasing = T),2])),10)
# }
# split = data$V2
# names(MC) = unique(split)
# 
# mat3 = mat2[order(match(rownames(mat2), data$V1)),] 
# pval3 = pval[order(match(rownames(pval), data$V1)),] 
# row.subsections <- c(length(cluster[[1]]),length(cluster[[2]]),length(cluster[[3]]),length(cluster[[4]]),
#                      length(cluster[[5]]),length(cluster[[6]]),length(cluster[[7]]),length(cluster[[8]]))
# row_split = data.frame(rep(c("MP_C1", "MP_C2", "MP_C3", "MP_C4", "MP_C5",
#                              "MP_C6", "MP_C7", "MP_C8"), row.subsections))
# pdf(paste0(path_result,"TCGA_MP_association_ImmuneResponse2.pdf"),width=25,height=19)
# set.seed(1129)
# Heatmap(mat3,name="Coefficient",
#         cluster_rows = F,
#         cluster_columns = T,
#         column_title = "Association between MP and Immune Response",
#         show_row_names = T,show_row_dend = T,
#         show_column_names = T,
#         column_names_rot=30,
#         cell_fun = function(j,i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         },
#         row_split = row_split,row_gap = unit(4, "mm"),
#         column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
#         column_title_gp = gpar(fontsize = 20,fontface="bold"),
#         right_annotation = rowAnnotation(textbox = anno_textbox(split, MC)))
# dev.off()
# 
# MC1 = list()
# MC1[[1]] = "Cold: Macrophage M2 Enriched -> High Hazard"
# MC1[[2]] = "Cold: Low Treg -> Low Hazard"
# MC1[[3]] = "Cold: Macrophage M2 Enriched, High ITH -> High Hazard"
# MC1[[4]] = "Cold: Treg Enriched -> High Hazard"
# MC1[[5]] = "Hot: Low M2, Low ITH -> Low Hazard"
# MC1[[6]] = "Hot: Treg Enriched -> High Hazard"
# MC1[[7]] = "Hot: Macrophage M1 Enriched,High LIL -> Low Hazard"
# MC1[[8]] = "Hot: Macrophage M2 Enriched -> High Hazard"
# names(MC1) = unique(split)
# 
# pdf(paste0(path_result,"TCGA_MP_association_ImmuneResponse3.pdf"),width=25,height=19)
# Heatmap(mat3,name="Coefficient",
#         cluster_rows = F,
#         cluster_columns = T,
#         column_title = "Association between MP and Immune Response",
#         show_row_names = T,show_row_dend = T,
#         show_column_names = T,
#         column_names_rot=30,
#         cell_fun = function(j,i, x, y, width, height, fill) {
#           grid.text(ifelse(pval[i,j]<0.05,"*",""), x, y, gp = gpar(fontsize = 10))
#         },
#         row_split = row_split,row_gap = unit(4, "mm"),
#         column_names_gp = grid::gpar(fontsize = 15,fontface="bold"),
#         column_title_gp = gpar(fontsize = 20,fontface="bold"),
#         right_annotation = rowAnnotation(textbox = anno_textbox(split, MC1),
#                                          annotation_name_gp= gpar(fontsize = 30,fontface="bold")),
#         )
# dev.off()
# 
# IM = read_excel(paste0(path_data,"immunomodulators.xlsx"),
#                         sheet = "Direct Relationship")[,1]
# IM = unlist(unname(as.vector(IM)))







  