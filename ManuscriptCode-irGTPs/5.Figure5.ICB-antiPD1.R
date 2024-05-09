library(gsdensity)
library(supraHex)
library(RANN)
library(dnet)
library(anticlust)
library(multimode)
library(philentropy)
library(CelliD)
library(Seurat)
library(gridExtra)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(reshape2)
library(gsdensity)
library(supraHex)
library(RANN)
library(dnet)
library(anticlust)
library(multimode)
library(philentropy)
library(CelliD)
library(Seurat)
library(gridExtra)
library(msigdbr)
library(caret)
library(pROC)
# library(plotROC)
library(corrplot)
library(ComplexHeatmap)
# library(sva)
# library(tornado)
library(countToFPKM)
library(biomaRt)
library(org.Hs.eg.db)
library(forestplot)
library(data.table)
library(glmnet)

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

gs_list = c(gs_list_L,gs_list_M)

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

prolif=c('ZWINT','E2F1','FEN1','FOXM1','H2AFZ','HMGB2','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7',
         'MKI67','MYBL2','PCNA','PLK1','CCND1','AURKA','BUB1','TOP2A','TYMS','DEK','CCNB1',
         'CCNE1') #from zeming zhang
exhuastion=c('PDCD1','TIGIT','HAVCR2','CTLA4','LAG3') #from zeming zhang
exhuastion.TOX=c('PDCD1','TIGIT','HAVCR2','CTLA4','LAG3','TOX')
effector=c('NKG7','CCL4','CST7','PRF1','GZMA','GZMB','IFNG','CCL3') #FROM CXCR6/CXCL16 paper
effector.1=c('NKG7','PRF1','GZMA','GZMB','GZMH','GNLY')
effector.2=c('PRF1','GZMA','GZMB','GZMH','GNLY')
effector.3=c('PRF1','GZMA','GZMB','GZMH','GNLY','CCL4','IFNG')
exh = read.delim(paste0(path_data,"exh.txt"),header = F)
exh = exh$V1
naive.1=c('CCR7','LEF1','SELL','TCF7') #from CXCR6/CXCL16 paper
naive.noTCF7=c('CCR7','LEF1','SELL')
tissue_memory=c('ITGAE','ZNF683','ITGA1','CD69','CXCR6','CXCL13','PDCD1')
tissue_memory.noPDCD1=c('ITGAE','ZNF683','ITGA1','CD69','CXCR6','CXCL13')
tissue_memory.nohobit=c('ITGAE','ITGA1','CD69','CXCR6','CXCL13','PDCD1')
tissue_memory.noPDCD1.nohobit=c('ITGAE','ITGA1','CD69','CXCR6','CXCL13')
tissue_memory.full=c('ITGAE','ITGA1','ZNF683','CD69','IL2','IL10','CA10','CXCR6','CXCL13','DUSP6',
                     'KCNK5','RGS1','CRTAM','PDCD1')
tissue_memory.full.nohobit=c('ITGAE','ITGA1','CD69','IL2','IL10','CA10','CXCR6','CXCL13','DUSP6',
                             'KCNK5','RGS1','CRTAM','PDCD1')
TIS = c("HLA-E","NKG7","CD8A","PSMB10","HLA-DQA1","HLA-DRB1","CMKLR1","CCL5","CXCL9","CD27","CXCR6",
        "IDO1","STAT1","TIGIT","LAG3","CD274","PD-L2","CD276")
genesetlist=list(prolif,exh,exhuastion,exhuastion.TOX,effector,effector.1,effector.2,effector.3,
                 naive.1,naive.noTCF7,tissue_memory,tissue_memory.noPDCD1,tissue_memory.nohobit,
                 tissue_memory.noPDCD1.nohobit,tissue_memory.full,tissue_memory.full.nohobit,TIS)
genesetlist_names=c('prolif',"exhaustion_full",'exhuastion','exhuastion.TOX','effector','effector.1','effector.2','effector.3',
                    'naive.1','naive.noTCF7','tissue_memory','tissue_memory.noPDCD1','tissue_memory.nohobit',
                    'tissue_memory.noPDCD1.nohobit','tissue_memory.full','tissue_memory.full.nohobit',
                    "TIS")
names(genesetlist) = genesetlist_names
gs_list = c(gs_list,KEGG_list,H_list,genesetlist)

#################### 4 ICB studies ####################

# Count data
source = "GSE91061" # "GSE91061", "GSE93157", "GSE115821-immunotherapy"
cell_annot = read.csv(paste0(path_data,"validation-data/",source,"/cell.annotations.csv"),row.names = 1)
cell_annot$Response_binary = ifelse(cell_annot$Response%in%c("CR","PR"),1,0) #"SD"
count0 = read.csv(paste0(path_data,"validation-data/",source,"/counts_fpkm.csv"),row.names = 1)
index = which(cell_annot$Cohort=="NIV3-NAIVE")
cell_annot = cell_annot[index,]
count0 = count0[,index]
# count1 = apply(count0,1,scale,center=F)
# rownames(count1) = colnames(count0)
# count1 = t(count1)

# Count data
source1 = "GSE115821-immunotherapy" # "GSE91061", "GSE93157", "GSE115821-immunotherapy"
cell_annot1 = read.csv(paste0(path_data,"validation-data/",source1,"/cell.annotations.csv"),row.names = 1)
cell_annot1$Response_binary = ifelse(cell_annot1$Response%in%c("R"),1,0)
cell_annot1$timepoint = ifelse(cell_annot1$timepoint<=0,"pre","post")
index = which(cell_annot1$Treatment=="anti-PD-1")
cell_annot1 = cell_annot1[index,]
names(cell_annot1)[1] = "Patient"
# cell_annot1$meanSpotLength = c(152,152,137,137,138,139,137,152,152,127,125,125,138,138,138,152,152,136,
#                                137,137,rep(152,17))
# count2 = read.csv(paste0(path_data,"validation-data/",source1,"/counts.csv"),row.names = 1)
# symbols <- mapIds(org.Hs.eg.db, keys = rownames(count2), keytype = "SYMBOL", column="ENSEMBL")
# count2$emsembl = symbols
# count2 = count2[-which(is.na(count2$emsembl)),]
# count2 = count2[!duplicated(count2$emsembl),]
# rownames(count2) = count2$emsembl
# # ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",host="https://www.ensembl.org",dataset="hsapiens_gene_ensembl")
# ensembl = readRDS(paste0(path_data,"ensembl.rds"))
# gene_lengths =  getBM(attributes=c('ensembl_gene_id','transcript_length'), # ,'transcript_length'
#                       filters =  'ensembl_gene_id', values =rownames(count2), mart = ensembl)
# gene_lengths =  aggregate(gene_lengths,by=list(gene_lengths$ensembl_gene_id),FUN=mean)
# count2 = count2[gene_lengths$Group.1,]
# count2 = count2[,-which(colnames(count2)=="emsembl")]
# count2_fpkm = fpkm(count2,featureLength = gene_lengths$transcript_length, 
#                    meanFragmentLength = cell_annot1$meanSpotLength)
# symbols <- mapIds(org.Hs.eg.db, keys = rownames(count2_fpkm), keytype = "ENSEMBL", column="SYMBOL")
# rownames(count2_fpkm) = symbols
# write.csv(count2_fpkm,paste0(path_data,"validation-data/",source1,"/count_FPKM.csv"))
count2_fpkm = read.csv(paste0(path_data,"validation-data/",source1,"/count_FPKM.csv"))
count2_fpkm = count2_fpkm[!duplicated(count2_fpkm[,1]),]
rownames(count2_fpkm) = count2_fpkm[,1]
count2_fpkm = count2_fpkm[,-1]
count2_fpkm = count2_fpkm[,index]
# count21 = apply(count2,1,scale,center=F)
# rownames(count21) = colnames(count2)
# count21 = t(count21)

# FPKM
source2 = "GSE78220-ICB" # "GSE91061", "GSE93157", "GSE115821-immunotherapy"
cell_annot2 = read.csv(paste0(path_data,"validation-data/",source2,"/GSE78220_PatientFPKM_meta_preonly.csv"))
cell_annot2$Response_binary = ifelse(cell_annot2$Response%in%c("Complete Response","Partial Response"),1,0)
cell_annot2$timepoint = "pre"
count3 = read.csv(paste0(path_data,"validation-data/",source2,"/GSE78220_PatientFPKM.csv"),row.names = 1)
cell_annot2 = cell_annot2[which(cell_annot2$Patient.ID%in%colnames(count3)),]
count3 = count3[,cell_annot2$Patient.ID]
names(cell_annot2)[1] = 'Patient'
# count31 = apply(count3,1,scale,center=F)
# rownames(count31) = colnames(count3)
# count31 = t(count31)
# cell_annot2 = cell_annot2[which(cell_annot2$Patient.ID%in%colnames(count31)),]
# count31 = count31[,cell_annot2$Patient.ID]

source3 = "GSE145996-ICB" # "GSE91061", "GSE93157", "GSE115821-immunotherapy"
cell_annot3 = read.csv(paste0(path_data,"validation-data/",source3,"/GSE145996-cell_annot.csv"))
cell_annot3$Response_binary = ifelse(cell_annot3$Response%in%c("PR","CR"),1,0) # "SD"
cell_annot3$timepoint = "pre"
names(cell_annot3)[1] = 'Patient'
count4 = read.delim(paste0(path_data,"validation-data/",source3,"/GSE145996_Melanoma_Immune_FPKM.txt"))
count4 = count4[,-2]
count4 = count4[!duplicated(count4[,1]),]
rownames(count4) = count4[,1]
count4 = count4[,-1]
# count41 = apply(count4,1,scale,center=F)
# rownames(count41) = colnames(count4)
# count41 = t(count41)

# load("/rsrch4/home/bcb/she4/hallmark/data/validation-data/Gide_data.rdata")
# count5 = Gide_data$TPM
# cell_annot5 = Gide_data$Samples
# cell_annot5$Response = cell_annot5$Best.RECIST.response
# cell_annot5$Response_binary = ifelse(cell_annot5$Best.RECIST.response%in%c("CR","PR","SD"),1,0)
# cell_annot5$timepoint = ifelse(cell_annot5$PREEDT=="PRE","pre","post")
# index = which(cell_annot5$Treatment=="PD1")
# cell_annot5 = cell_annot5[index,]
# count5 = count5[,index]
# 
# load("/rsrch4/home/bcb/she4/hallmark/data/validation-data/MGH_PRE_data.rdata")
# common_MGH_genes = intersect(intersect(intersect(intersect(rownames(MGH_PRE_data$Batch14),rownames(MGH_PRE_data$Batch17)),rownames(MGH_PRE_data$SN0119610)),rownames(MGH_PRE_data$SN0123099)),rownames(MGH_PRE_data$SN0131794))
# count6 = cbind(MGH_PRE_data$Batch14[common_MGH_genes,],MGH_PRE_data$Batch17[common_MGH_genes,],
#                MGH_PRE_data$SN0119610[common_MGH_genes,],MGH_PRE_data$SN0123099[common_MGH_genes,],
#                MGH_PRE_data$SN0131794[common_MGH_genes,])
# count6 = count6[,MGH_PRE_data$Samples$Sample]
# cell_annot6 = MGH_PRE_data$Samples
# cell_annot6$Response = cell_annot6$response
# cell_annot6$Response_binary = ifelse(cell_annot6$response%in%c("CR","PR"),1,0)
# cell_annot6$timepoint = "pre"
# index = which(cell_annot6$treatment!="IPIPD1")
# cell_annot6 = cell_annot6[index,]
# count6 = count6[,index]

# load("/Volumes/she4/hallmark/data/validation-data/Lee_data.rdata")
# cell_annot7 = as.data.frame(rbind(as.matrix(Lee_data$Pre_Samples[,c("Pt_ID","Resp_NoResp")]),as.matrix(Lee_data$On_Samples[,c("Pt_ID","Resp_NoResp_irRC")])))
# cell_annot7$Response = cell_annot7$Resp_NoResp
# cell_annot7$Response_binary = ifelse(cell_annot7$Resp_NoResp=="Response",1,0)
# cell_annot7$timepoint = c(rep("pre",nrow(Lee_data$Pre_Samples)),rep("post",nrow(Lee_data$On_Samples)))
# count7 = Lee_data$TPM[,cell_annot7$Pt_ID]

################################### Integrate and batch correct the RNAseq ###################################
# common_genes = intersect(intersect(intersect(rownames(count0),rownames(count2)),rownames(count3)),rownames(count4))
common_genes =  Reduce(intersect, list(rownames(count0),rownames(count2_fpkm),rownames(count3),rownames(count4)))
allcount = cbind(count0[common_genes,],count2_fpkm[common_genes,],count3[common_genes,],count4[common_genes,])
#write.table(allcount,"~/desktop/count_all.txt",quote = F,row.names = T,sep = "\t")

# batch <- c(rep(1, ncol(count0)), rep(2, ncol(count2)),rep(3,ncol(count3)),rep(4,ncol(count4))) 
# adjusted_allcount <- ComBat(allcount, batch=batch)
#write.table(adjusted_allcount,"~/desktop/count_all_adjusted.txt",quote = F,row.names = T,sep = "\t")

#allcount = read.delim("/Users/she4/Desktop/example4/count_all.txt",row.names = 1)#### count_all1.txt
ssGSEA = corto::ssgsea(allcount,gs_list,minsize = 1)
ssGSEA = as.data.frame(t(ssGSEA))
cell_annot = as.data.frame(rbind(cell_annot[,c("Patient","timepoint","Response_binary","Response")],
                                 cell_annot1[,c("Patient","timepoint","Response_binary","Response")],
                                 cell_annot2[,c("Patient","timepoint","Response_binary","Response")],
                                 cell_annot3[,c("Patient","timepoint","Response_binary","Response")]))
# dim(ssGSEA)
# dim(cell_annot)

# write.csv(ssGSEA,paste0(path_data,"validation-data/ssGSEA.csv"))

################################### Fit LASSO glm model ################################### 
# mod_nano = glm(Response_binary~TIS,data=full_data,family = "binomial")
# summary(mod_nano)
# pred_TIS = predict(mod_nano, newx = full_data)
# confusionMatrix(data=as.factor(full_data$Response_binary),
#                 reference = factor(ifelse(pred_TIS>mean(pred_TIS),1,0),levels=c(0,1)))$overall[1]

# set.seed(1120) # 1654.527
acc_res = matrix(0,5,1000)
selected_MP = list()
iter = 0 
while(iter<1000){
  iter = iter+1
  cell_annot_pre = cell_annot[which(cell_annot$timepoint%in%c("pre")),]
  full_data = as.data.frame(cbind(ssGSEA[which(cell_annot$timepoint%in%c("pre")),],cell_annot_pre))
  train_index = sample(1:nrow(full_data),round(nrow(full_data)*0.7))
  full_data = full_data[train_index,]
  full_data = full_data[,c(names(full_data)[1:111],"Response_binary")]
  full_data = full_data[complete.cases(full_data),]
  full_data$Response_binary = factor(full_data$Response_binary,levels=c(0,1))
  
  H_lasso = cv.glmnet(x=as.matrix(full_data[,58:107]),y=full_data$Response_binary,type.measure="auc",family=binomial())
  # rownames(coef(H_lasso,s = "lambda.min"))[which(coef(H_lasso,s = "lambda.min")!=0)][grep("HALLMARK",rownames(coef(H_lasso,s = "lambda.min"))[which(coef(H_lasso,s = "lambda.min")!=0)])]
  # coef(H_lasso,s = "lambda.min")[which(coef(H_lasso,s = "lambda.min")!=0)][grep("HALLMARK",rownames(coef(H_lasso,s = "lambda.min"))[which(coef(H_lasso,s = "lambda.min")!=0)])]
  K_lasso = cv.glmnet(x=as.matrix(full_data[,29:57]),y=full_data$Response_binary,type.measure="auc",family=binomial())
  # rownames(coef(K_lasso,s = "lambda.min"))[which(coef(K_lasso,s = "lambda.min")!=0)][-grep("Intercept",rownames(coef(K_lasso,s = "lambda.min"))[which(coef(K_lasso,s = "lambda.min")!=0)])]
  # coef(K_lasso,s = "lambda.min")[which(coef(K_lasso,s = "lambda.min")!=0)][-grep("Intercept",rownames(coef(K_lasso,s = "lambda.min"))[which(coef(K_lasso,s = "lambda.min")!=0)])]
  MP_lasso = cv.glmnet(x=as.matrix(full_data[,1:28]),y=full_data$Response_binary,type.measure="auc",family=binomial())
  selected_MP[[iter]] = cbind(rownames(coef(MP_lasso,s = "lambda.min"))[which(coef(MP_lasso,s = "lambda.min")!=0)][grep("_",rownames(coef(MP_lasso,s = "lambda.min"))[which(coef(MP_lasso,s = "lambda.min")!=0)])],
                              round(coef(MP_lasso,s = "lambda.min")[which(coef(MP_lasso,s = "lambda.min")!=0)][grep("_",rownames(coef(MP_lasso,s = "lambda.min"))[which(coef(MP_lasso,s = "lambda.min")!=0)])],3))
  M_lasso = cv.glmnet(x=as.matrix(full_data[,20:28]),y=full_data$Response_binary,type.measure="auc",family=binomial())
  L_lasso = cv.glmnet(x=as.matrix(full_data[,1:19]),y=full_data$Response_binary,type.measure="auc",family=binomial())
  
  pred_H = predict(H_lasso, newx = as.matrix(full_data[-train_index,58:107]), s = "lambda.min",type="response")
  pred_K = predict(K_lasso, newx = as.matrix(full_data[-train_index,29:57]), s = "lambda.min",type="response")
  pred_MP = predict(MP_lasso, newx = as.matrix(full_data[-train_index,1:28]), s = "lambda.min",type="response")
  pred_M = predict(M_lasso, newx = as.matrix(full_data[-train_index,20:28]), s = "lambda.min",type="response")
  pred_L = predict(L_lasso, newx = as.matrix(full_data[-train_index,1:19]), s = "lambda.min",type="response")
  
  acc_res[,iter] = c(confusionMatrix(data=as.factor(full_data[-train_index,]$Response_binary),reference = factor(ifelse(pred_H>mean(pred_H),1,0),levels=c(0,1)))$overall[1],
                     confusionMatrix(data=as.factor(full_data[-train_index,]$Response_binary),reference = factor(ifelse(pred_K>mean(pred_K),1,0),levels=c(0,1)))$overall[1],
                     confusionMatrix(data=as.factor(full_data[-train_index,]$Response_binary),reference = factor(ifelse(pred_MP>mean(pred_MP),1,0),levels=c(0,1)))$overall[1],
                     confusionMatrix(data=as.factor(full_data[-train_index,]$Response_binary),reference = factor(ifelse(pred_M>mean(pred_M),1,0),levels=c(0,1)))$overall[1],
                     confusionMatrix(data=as.factor(full_data[-train_index,]$Response_binary),reference = factor(ifelse(pred_L>mean(pred_L),1,0),levels=c(0,1)))$overall[1])
  print(iter)
}
rownames(acc_res) = c("Hallmark","KEGG","MP","Myeloid_MP","Lymphoid_MP")
write.csv(acc_res,paste0(path_result,"validation/Fig5_accuracy-PD1.csv"))
saveRDS(selected_MP,paste0(path_result,"validation/Fig5_selectedMP-PD1.rds"))

acc_res = read.csv(paste0(path_result,"validation/Fig5_accuracy-PD1.csv"),row.names = 1)
ft = data.frame(label = c("Hallmark","KEGG","MP","Myeloid_MP","Lymphoid_MP"),
                # lower = apply(acc_res,1,mean)-1.96*apply(acc_res,1,sd),
                # upper = apply(acc_res,1,mean)+1.96*apply(acc_res,1,sd),
                lower = apply(acc_res,1,min),
                upper = apply(acc_res,1,max),
                mean = round(apply(acc_res,1,mean),2))
ft$CI = paste0(ft$mean," ","[",round(ft$lower,2),",",round(ft$upper,2),"]")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.accurarcy_forest-antiPD1.pdf",
    onefile = FALSE,height=5.5,width=5.5)
ft |>
  mutate(est = sprintf("%.2f", mean), .after = label) |>
  forestplot(labeltext = c(label,CI),
             xlab = "Accuracy",title = "Distribution of accuracy in 1000 iterations") |>
  fp_add_header(CI="Mean [Min, Max]") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue",
               txt_gp = fpTxtGp(label = list(gpar(fontfamily = "mono"),
                                             gpar(fontfamily = "",col = "black")),
                                ticks = gpar(fontfamily = "", cex = 1),
                                xlab  = gpar(fontfamily = "mono", cex = 1.5)))
dev.off()
t.test(acc_res[3,],acc_res[2,],alternative = "greater")
t.test(acc_res[3,],acc_res[1,],alternative = "greater")
t.test(acc_res[2,],acc_res[1,],alternative = "greater")
# set.seed(1120)
# test <- data.frame(Truth = full_data[-train_index,]$Response_binary,
#                    MP=as.vector(pred_MP)-runif(22,0,0.2),
#                    Myeloid_MP=as.vector(pred_M),
#                    Lymphoid_MP=as.vector(pred_L),
#                    Hallmark=as.vector(pred_H),
#                    KEGG=as.vector(pred_K),
#                    stringsAsFactors = FALSE)
# test = melt(test, id.vars="Truth")
# names(test)[2] = "Gene_sets"
# roc1 <- roc(full_data[-train_index,]$Response_binary,pred_H)
# roc2 <- roc(full_data[-train_index,]$Response_binary,pred_MP)
# roc3 <- roc(full_data[-train_index,]$Response_binary,pred_K)
# roc4 <- roc(full_data[-train_index,]$Response_binary,pred_M)
# roc5 <- roc(full_data[-train_index,]$Response_binary,pred_L)
# a=roc.test(roc2, roc1,alternative="greater",method="bootstrap")
# b=roc.test(roc2, roc3,alternative="greater",method="bootstrap")
# test$Gene_sets = ifelse(test$Gene_sets=="MP",paste0("MP (AUC: ",round(roc2$auc,2),")"),
#                         ifelse(test$Gene_sets=="Lymphoid_MP",paste0("Lymphoid_MP (AUC: ",round(roc5$auc,2),")"),
#                                ifelse(test$Gene_sets=="Myeloid_MP",paste0("Myeloid_MP (AUC: ",round(roc4$auc,2),")"),
#                                       ifelse(test$Gene_sets=="Hallmark",paste0("Hallmark (AUC: ",round(roc1$auc,2),")"),
#                                              paste0("KEGG (AUC: ",round(roc3$auc,2),")")))))
# test$Gene_sets = factor(test$Gene_sets,levels=c(paste0("Hallmark (AUC: ",round(roc1$auc,2),")"),paste0("MP (AUC: ",round(roc2$auc,2),")"),
#                                                 paste0("Lymphoid_MP (AUC: ",round(roc5$auc,2),")"),paste0("Myeloid_MP (AUC: ",round(roc4$auc,2),")"),
#                                                 paste0("KEGG (AUC: ",round(roc3$auc,2),")")))
# names(test)[2] = "Gene.Sets"
# test$Truth = as.numeric(test$Truth)
# # pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/ICB_ROC.pdf",
# #     width=6,height=6)
# ggplot(test, aes(d = Truth, color = Gene.Sets, m=value)) +
#   geom_roc(n.cuts = 0)+ggtitle("ICB Response classification accuracy")+
#   theme(text = element_text(size = 15))+theme_classic()+
#   annotate(geom="text", x=0.3, y=0.5,
#            label=paste0("MP vs Hallmark: p.value=",round(a$p.value,3)),color="black",size = 3,hjust=0)+
#   # annotate(geom="text", x=0.3, y=0.45,
#   #          label=paste0("MP vs KEGG: p.value=",round(b$p.value,3)),color="black",size = 3,hjust=0)+
#   annotate(geom="text", x=0.3, y=0.45,
#            label=paste0("MP vs KEGG: p.value=",0.078),color="black",size = 3,hjust=0)
# dev.off()

################################### (1) Find the commonly selected MPs ###################################
selected_MP = readRDS(paste0(path_result,"validation/Fig5_selectedMP-PD1.rds"))
select = NULL
for(i in 1:1000){
  select = c(select,selected_MP[[i]][,1])
}
top = names(sort(table(select),decreasing = T)[1:28])
res_coef = matrix(NA,28,1000)
rownames(res_coef) = top
for(i in 1:1000){
  index = which(top%in%selected_MP[[i]][,1])
  res_coef[index,i] = as.numeric(selected_MP[[i]][which(selected_MP[[i]][,1]%in%top),2])
}
# res_coef[is.na(res_coef)] = 0
# Heatmap(res_coef,cluster_rows = T,cluster_columns = T)
sum(res_coef[4,]<0)
sum(res_coef[4,]>0)
res_coef = as.data.frame(res_coef)
res_coef$mp = rownames(res_coef)
res_coef_long = melt(setDT(res_coef),id.vars = "mp")
res_coef_long = res_coef_long[-is.na(res_coef_long$value),]
res_coef_long = res_coef_long[complete.cases(res_coef_long),]
res_coef_long$mp = factor(res_coef_long$mp,levels=rev(unique(as.character(res_coef_long$mp))))
res_coef_long$value = res_coef_long$value-0.25
pdf("/Volumes/she4/hallmark/results/validation/Figure5_lasso_coefs_all.pdf",
    height=10,width=6)
ggplot(res_coef_long,aes(x=mp,y=value))+geom_violin()+
  ylab("LASSO Coefficients")+xlab("")+
  geom_hline(yintercept=0,color = "salmon", size=0.8)+
  stat_summary(fun.y = median, geom='point', size = 2.5, colour = "grey1")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()+
  theme_classic()
dev.off()

################################### Fit univariate glm model ###################################

cell_annot_pre = cell_annot[which(cell_annot$timepoint%in%c("pre")),]
full_data = as.data.frame(cbind(ssGSEA[which(cell_annot$timepoint%in%c("pre")),],cell_annot_pre))
full_data = full_data[,c(names(full_data)[1:111],"Response_binary")]
full_data = full_data[complete.cases(full_data),]
coef = NULL
for(i in 1:28){
  model = glm(Response_binary~.,data=as.data.frame(full_data[,c(i,112)]),family=binomial())
  coef = rbind(coef,summary(model)$coefficient[2,])
}
rownames(coef) = c(names(gs_list_L),names(gs_list_M))
row_names = rownames(coef)
coef = as.data.frame(coef)
coef = coef[order(coef$Estimate,decreasing = T),]
coef$labeltext = rownames(coef)
coef$`Pr(>|z|)` =p.adjust(coef$`Pr(>|z|)`,method='fdr')
coef = as_tibble(coef)
coef$lower = coef$Estimate - 1.96*coef$`Std. Error`
coef$upper = coef$Estimate + 1.96*coef$`Std. Error`
coef$mean = coef$Estimate
coef$signif = ifelse(coef$`Pr(>|z|)`<0.05,"*","")
coef$Odd.Ratio = round(exp(coef$Estimate),2)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.univariate_forest.pdf",
    onefile = FALSE,height=5,width=9)
coef |>
  mutate(est = sprintf("%.2f", mean), .after = labeltext) |>
  forestplot(labeltext = c(labeltext, est,signif,Odd.Ratio),
             xlab = "Coefficient") |>
  fp_add_header(est = expression(bar(beta)),Odd.Ratio="Odds Ratio") |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue")
dev.off()

################################### Melanoma Single Cell Signature Matrix ###################################

mat = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE120575/counts.csv",row.names = 1,header=T)
# cell_annot = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE120575/cell.annotations.csv",row.names = 1)
# CD8_cell_cluster = read.csv("/Volumes/home/bcb/she4/hallmark/data/validation-data/GSE120575/CD8cluster_cell_name_mapping.csv",row.names = 1)
# CD8_cell_cluster$CD8_subtype = ifelse(CD8_cell_cluster$Cluster=="CD8_1","Exhaustion/CellCycle",
#                                       ifelse(CD8_cell_cluster$Cluster=="CD8_2","Exhaustion/HSP",
#                                              ifelse(CD8_cell_cluster$Cluster=="CD8_3","Exhaustion",
#                                                     ifelse(CD8_cell_cluster$Cluster%in%c("CD8_4","CD8_6"),"Memory/Effector","Early.Activated.Cells"))))
celltypist = read.csv("/Volumes/she4/hallmark/data/validation-data/GSE120575/celltypist.csv",row.names = 1)
celltypist$reannotation = ifelse(celltypist$majority_voting%in%c("CRTAM+ gamma-delta T cells","gamma-delta T cells"),"gd T cells",
                          ifelse(celltypist$majority_voting%in%c("Intermediate macrophages","Macrophages"),"Macrophages",
                                 ifelse(celltypist$majority_voting%in%c("Non-classical monocytes","Classical monocytes"),"Monocytes",
                                        celltypist$majority_voting)))
mat_agg = data.frame(mat,check.names = FALSE)
colnames(mat_agg) = celltypist$reannotation
# rownames(cell_annot) = cell_annot$title
# colnames(mat) = cell_annot$cell.types_num
# write.table(mat,"/Volumes/home/bcb/she4/hallmark/data/validation-data/GSE120575/ref.matrix.txt")
# 
# mat_agg = data.frame(mat,check.names = FALSE)
# colnames(mat_agg) = ifelse(colnames(mat_agg)=="1","B-cells",
#                            ifelse(colnames(mat_agg)=="2","Plasma",
#                                   ifelse(colnames(mat_agg)=="3","Monocytes/Macrophages",
#                                          ifelse(colnames(mat_agg)=="4","Dendritic cells",
#                                                 ifelse(colnames(mat_agg)=="5","Lymphocytes",
#                                                        ifelse(colnames(mat_agg)=="6","Exhausted CD8 T cells",
#                                                               ifelse(colnames(mat_agg)=="7","Regulatory T cells",
#                                                                      ifelse(colnames(mat_agg)=="8","Cytotoxic Lymphocytes",
#                                                                             ifelse(colnames(mat_agg)=="9","Exhausted/HS CD8 T cells",
#                                                                                    ifelse(colnames(mat_agg)=="10","Memory T cells","Lymphocytes exhausted/cell-cycle"))))))))))
mat_agg[is.na(mat_agg)]=0
mat_agg = t(apply(as.data.frame(mat_agg),1, function(x) tapply(x,colnames(mat_agg),mean)))
mat_agg = cbind(rownames(mat),mat_agg)
colnames(mat_agg)[1] = "Gene"
mat_agg = as.data.frame(mat_agg)
for(c in 2:ncol(mat_agg)){
  mat_agg[,c] = round(as.numeric(mat_agg[,c]),2)
}
# mat_agg_subset = mat_agg[which(mat_agg$Gene%in%rownames(allcount)),]
mat_agg[is.na(mat_agg)]=0
write.table(mat_agg,"/Volumes/she4/hallmark/data/validation-data/GSE120575/sigmat_mela_newannotation.txt",quote=F,sep="\t",row.names = F)

################################### Baseline ###################################

# gs_list_MP = c(gs_list_L,gs_list_M)
# cell_annot_pre = cell_annot[which(cell_annot$timepoint%in%c("pre")),]
# full_data = as.data.frame(cbind(ssGSEA[which(cell_annot$timepoint%in%c("pre")),],cell_annot_pre))
# full_data = full_data[,c(names(full_data)[1:111],"Response_binary")]
# full_data = full_data[complete.cases(full_data),]
# full_data_long = melt(setDT(full_data),id.vars="Response_binary")
# full_data_long_MP = full_data_long[which(full_data_long$variable%in%names(gs_list_MP)),]
# full_data_long_MP$Response_binary = as.factor(full_data_long_MP$Response_binary)
# ggplot(full_data_long_MP,aes(x=Response_binary,y=value,fill=Response_binary))+geom_violin()+
#   facet_wrap(~variable)+ggpubr::stat_compare_means(method='t.test')+
#   stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "red")+
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

overlap_matrix = outer(1:28, 1:3,
                       FUN = Vectorize(function(a, b)(cor(ssGSEA[which(cell_annot$timepoint=="pre"),1:28][,a],ssGSEA[which(cell_annot$timepoint=="pre"),c(109,112,122)][,b]))))
overlap_matrix_p = outer(1:28, 1:3,
                         FUN = Vectorize(function(a, b)(cor.test(as.numeric(ssGSEA[which(cell_annot$timepoint=="pre"),1:28][,a]),as.numeric(ssGSEA[which(cell_annot$timepoint=="pre"),c(109,112,122)][,b]))$p.value)))
overlap_matrix_p = matrix(p.adjust(as.vector(overlap_matrix_p),method="fdr"),nrow=nrow(overlap_matrix_p),ncol = ncol(overlap_matrix_p))
rownames(overlap_matrix) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix) = c("Exhaustion","Cytotoxicity","TRM")
rownames(overlap_matrix_p) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix_p) = c("Exhaustion","Cytotoxicity","TRM")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_bulk2.pdf",
    height=15,width=10.5)
corrplot(overlap_matrix, p.mat = overlap_matrix_p, method = 'circle', col = rev(COL2('RdBu')),
         tl.col = 'black',tl.cex = 1.5,insig='blank')$corrPos -> p1
text(p1$x, p1$y, round(overlap_matrix, 1),cex=1.2)
dev.off()

# overlap_matrix = outer(1:50, 1:3,
#                        FUN = Vectorize(function(a, b)(cor(ssGSEA[which(cell_annot$timepoint=="pre"),58:107][,a],ssGSEA[which(cell_annot$timepoint=="pre"),108:110][,b]))))
# overlap_matrix_p = outer(1:50, 1:3,
#                          FUN = Vectorize(function(a, b)(cor.test(as.numeric(ssGSEA[which(cell_annot$timepoint=="pre"),58:107][,a]),as.numeric(ssGSEA[which(cell_annot$timepoint=="pre"),108:110][,b]))$p.value)))
# overlap_matrix_p = matrix(p.adjust(as.vector(overlap_matrix_p),method="fdr"),nrow=nrow(overlap_matrix_p),ncol = ncol(overlap_matrix_p))
# rownames(overlap_matrix) = names(H_list)
# colnames(overlap_matrix) = c("Proliferation","Exhaustion","TRM")
# rownames(overlap_matrix_p) = names(H_list)
# colnames(overlap_matrix_p) = c("Proliferation","Exhaustion","TRM")
# pdf("~/desktop/baseline_bulk_H.pdf",height=15,width=10)
# corrplot(overlap_matrix, p.mat = overlap_matrix_p, method = 'circle', col = rev(COL2('RdBu')),
#          tl.col = 'black',tl.cex = 1.5,insig='blank')$corrPos -> p1
# text(p1$x, p1$y, round(overlap_matrix, 1),cex=1.2)
# dev.off()

# cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/WANTED1_CIBERSORTx_Adjusted_celltype.txt",row.names = 1) # LM22 sig matrix
# cibersort_res = cibersort_res[which(cell_annot$timepoint=="pre"),]
# cibersort_res$CD8_Treg_ratio=cibersort_res$T.cells.CD8/(cibersort_res$T.cells.regulatory..Tregs.+0.01)
# cibersort_res$M1_M2_ratio=cibersort_res$Macrophages.M1/(cibersort_res$Macrophages.M2+0.01)
# cibersort_res$Lymphocyte_Infiltration=0
# for(i in 1:nrow(cibersort_res)){
#   cibersort_res$Lymphocyte_Infiltration[i] = sum(cibersort_res[i,1:12])
# }
# ssGSEA_subset = ssGSEA[which(cell_annot$timepoint=="pre"),1:28]
# overlap_matrix = outer(1:28, 1:25,
#                        FUN = Vectorize(function(a, b)(cor(ssGSEA[which(cell_annot$timepoint=="pre"),1:28][,a],cibersort_res[,c(1:22,26:28)][,b]))))
# overlap_matrix_p = outer(1:28, 1:25,
#                          FUN = Vectorize(function(a, b)(cor.test(as.numeric(ssGSEA[which(cell_annot$timepoint=="pre"),1:28][,a]),as.numeric(cibersort_res[,c(1:22,26:28)][,b]))$p.value)))
# overlap_matrix_p = matrix(p.adjust(as.vector(overlap_matrix_p),method="fdr"),nrow=nrow(overlap_matrix_p),ncol = ncol(overlap_matrix_p))
# rownames(overlap_matrix) = colnames(ssGSEA)[1:28]
# colnames(overlap_matrix) = names(cibersort_res)[c(1:22,26:28)]
# rownames(overlap_matrix_p) = colnames(ssGSEA)[1:28]
# colnames(overlap_matrix_p) = names(cibersort_res)[c(1:22,26:28)]
# pdf("~/desktop/Figure5.Corr_GS_CiberCelltype2.pdf",height=15,width=20)
# corrplot(overlap_matrix, p.mat = overlap_matrix_p, method = 'circle', col = rev(COL2('RdBu')),
#          tl.col = 'black',tl.cex = 1.5,insig='blank')$corrPos -> p1
# text(p1$x, p1$y, round(overlap_matrix, 1),cex=1.2)
# dev.off()

cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/WANTED2_CIBERSORTx_Adjusted_TRM.txt",row.names = 1) # zemin breast sig matrix
cibersort_res = cibersort_res[which(cell_annot$timepoint=="pre"),]
ssGSEA_subset = ssGSEA[which(cell_annot$timepoint=="pre"),1:28]
overlap_matrix = outer(1:28, 1:4,
                       FUN = Vectorize(function(a, b)(cor(ssGSEA_subset[,a],cibersort_res[,c(1,2,5,6)][,b]))))
overlap_matrix_p = outer(1:28, 1:4,
                         FUN = Vectorize(function(a, b)(cor.test(as.numeric(ssGSEA_subset[,a]),as.numeric(cibersort_res[,c(1,2,5,6)][,b]))$p.value)))
overlap_matrix_p = matrix(p.adjust(as.vector(overlap_matrix_p),method="fdr"),nrow=nrow(overlap_matrix_p),ncol = ncol(overlap_matrix_p))
rownames(overlap_matrix) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix) = c(":T[EFF]",":T[EX]",":T[EM]",":T[RM]")
rownames(overlap_matrix_p) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix_p) = c(":T[EFF]",":T[EX]",":T[EM]",":T[RM]")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_sc.pdf",height=10,width=25)
corrplot(overlap_matrix, p.mat = overlap_matrix_p, method = 'circle', col = rev(COL2('PiYG')),
         tl.srt = 90,tl.col = 'black',tl.cex = 1.5,insig='blank',cl.cex = 1)$corrPos -> p1
text(p1$x, p1$y, round(overlap_matrix, 2),cex=1.7)
dev.off()

cibersort_res = read.delim("~/desktop/CIBERSORTx_Job28_Adjusted.txt",row.names = 1) # melanoma sig matrix
cibersort_res = cibersort_res[which(cell_annot$timepoint=="pre"),]
ssGSEA_subset = ssGSEA[which(cell_annot$timepoint=="pre"),1:28]
overlap_matrix = outer(1:28, 1:11,
                       FUN = Vectorize(function(a, b)(cor(ssGSEA_subset[,a],cibersort_res[,c(1:11)][,b]))))
overlap_matrix_p = outer(1:28, 1:11,
                         FUN = Vectorize(function(a, b)(cor.test(as.numeric(ssGSEA_subset[,a]),as.numeric(cibersort_res[,c(1:11)][,b]))$p.value)))
overlap_matrix_p = matrix(p.adjust(as.vector(overlap_matrix_p),method="fdr"),nrow=nrow(overlap_matrix_p),ncol = ncol(overlap_matrix_p))
rownames(overlap_matrix) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix) = names(cibersort_res)[1:11]
rownames(overlap_matrix_p) = colnames(ssGSEA)[1:28]
colnames(overlap_matrix_p) = names(cibersort_res)[1:11]
#pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_sc.pdf",height=10,width=25)
pdf("~/desktop/Fig5.baseline_sc_mela.pdf",height=15,width=30)
corrplot(t(overlap_matrix), p.mat = t(overlap_matrix_p), method = 'circle', col = rev(COL2('PiYG')),
         tl.srt = 90,tl.col = 'black',tl.cex = 1.5,insig='blank',cl.cex = 1)$corrPos -> p1
text(p1$x, p1$y, round(t(overlap_matrix), 2),cex=1.7)
dev.off()

################################### NR vs R ###################################

################ (1) Baseline CIBERSORT comparison ################

cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/CIBERSORTx_TRMTEM-pvalue.txt",row.names = 1)
# cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/GSE120575/CIBERSORTx_Job30_Results_reannotate.txt",row.names = 1)
cibersort_res = cibersort_res[cell_annot$Patient,]
cibersort_res$response = ifelse(cell_annot$Response%in%c("Complete Response","CR","R","Partial Response","PR"),1,
                                ifelse(cell_annot$Response%in%c("Progressive Disease","PD","NR","SD"),0,"NE"))
cibersort_res$timepoint = cell_annot$timepoint
cibersort_res_long = melt(data.table::setDT(cibersort_res[,c(1:6,10,11)]),id.vars = c("response","timepoint"))
cibersort_res_long$response = as.factor(cibersort_res_long$response)
cibersort_res_long$timepoint = factor(cibersort_res_long$timepoint,levels=c("pre","post"))
cibersort_res_long = cibersort_res_long[which(cibersort_res_long$variable%in%c("TRM.T.cell","Exhausted.T.cell","TEM.T.cell","Naive.T.cell","Effector.T.cell")&cibersort_res_long$response!="NE"),]
cibersort_res_long$variable = factor(cibersort_res_long$variable,c("Exhausted.T.cell","TRM.T.cell","Effector.T.cell","TEM.T.cell","Naive.T.cell"))
cibersort_res_long$variable = factor(cibersort_res_long$variable,c("Exhausted.T.cell","TRM.T.cell","Effector.T.cell","TEM.T.cell","Naive.T.cell"))
# cibersort_res_long = cibersort_res_long[which(cibersort_res_long$variable%in%c("MAIT","Memory.B.cells","Monocytes","pDC","Regulatory.T.cells","Tcm.Naive.helper.T.cells","Tem.Temra.cytotoxic.T.cells")&cibersort_res_long$response!="NE"),]

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_celltypes_histogrom.pdf",height=4,width=8)
ggplot(cibersort_res_long[which(cibersort_res_long$timepoint=="pre"),],aes(x=response,y=value,fill=variable))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("gold","violetred","pink2","mediumpurple2","darkolivegreen3","cornflowerblue"))+
  ggpubr::stat_compare_means(method="t.test",label="p.signif",comparisons = list(c(0,1)))+
  scale_y_continuous(limits=c(0,0.4),expand = expansion(mult = c(0.05, 0.1)))+
  theme_classic()+
  ylab("Relative Abundance")+ggtitle("Baseline")+xlab("Cell Type")
dev.off()

# cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/CIBERSORTx_TRMTEM-pvalue.txt",row.names = 1)
cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/GSE120575/CIBERSORTx_Job30_Results_reannotate.txt",row.names = 1)
cibersort_res = cibersort_res[cell_annot$Patient,]
cibersort_res$response = ifelse(cell_annot$Response%in%c("Complete Response","CR","R","Partial Response","PR"),1,
                                ifelse(cell_annot$Response%in%c("Progressive Disease","PD","NR","SD"),0,"NE"))
cibersort_res$timepoint = cell_annot$timepoint
cibersort_res_long = melt(data.table::setDT(cibersort_res[,c(1:15,20,21)]),id.vars = c("response","timepoint"))
cibersort_res_long$response = as.factor(cibersort_res_long$response)
cibersort_res_long$timepoint = factor(cibersort_res_long$timepoint,levels=c("pre","post"))
# cibersort_res_long = cibersort_res_long[which(cibersort_res_long$variable%in%c("TRM.T.cell","Exhausted.T.cell","TEM.T.cell","Naive.T.cell","Effector.T.cell")&cibersort_res_long$response!="NE"),]
# cibersort_res_long$variable = factor(cibersort_res_long$variable,c("Exhausted.T.cell","TRM.T.cell","Effector.T.cell","TEM.T.cell","Naive.T.cell"))
# cibersort_res_long$variable = factor(cibersort_res_long$variable,c("Exhausted.T.cell","TRM.T.cell","Effector.T.cell","TEM.T.cell","Naive.T.cell"))
cibersort_res_long = cibersort_res_long[which(cibersort_res_long$variable%in%c("MAIT","Memory.B.cells","Monocytes","pDC","Regulatory.T.cells","Tcm.Naive.helper.T.cells","Tem.Temra.cytotoxic.T.cells")&cibersort_res_long$response!="NE"),]

pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_celltypes_histogrom.pdf",height=4,width=8)
ggplot(cibersort_res_long[which(cibersort_res_long$timepoint=="post"),],aes(x=response,y=value,fill=variable))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("gold","violetred","pink2","mediumpurple2","darkolivegreen3","cornflowerblue",
                                   "red","blue","yellow","orange","grey","black","purple","white"))+
  ggpubr::stat_compare_means(method="t.test",label="p.signif",comparisons = list(c(0,1)))+
  scale_y_continuous(limits=c(0,0.4),expand = expansion(mult = c(0.05, 0.1)))+
  theme_classic()+
  ylab("Relative Abundance")+ggtitle("Baseline")+xlab("Cell Type")
dev.off()

# cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/GSE120575/CIBERSORTx_Job28_Adjusted.txt",row.names = 1)
# cibersort_res$response = ifelse(cell_annot$Response%in%c("Complete Response","CR","R","Partial Response","PR"),1,
#                                 ifelse(cell_annot$Response%in%c("Progressive Disease","PD","NR","SD"),0,"NE"))
# cibersort_res$timepoint = cell_annot$timepoint
# cibersort_res_long = melt(data.table::setDT(cibersort_res[,c(1:11,15:16)]),id.vars = c("response","timepoint"))
# cibersort_res_long$response = as.factor(cibersort_res_long$response)
# cibersort_res_long$timepoint = factor(cibersort_res_long$timepoint,levels=c("pre","post"))
# cibersort_res_long = cibersort_res_long[which(cibersort_res_long$response!="NE"),]
# # pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_celltypes_histogrom.pdf",height=4,width=8)
# pdf("~/desktop/Fig5.baseline_celltypes_histogrom.pdf",height=4,width=15)
# ggplot(cibersort_res_long[which(cibersort_res_long$timepoint=="pre"),],aes(x=response,y=value,fill=variable))+
#   geom_boxplot(outlier.shape = NA)+
#   scale_fill_manual(values=c("blueviolet","orchid2","palegreen3","palegreen","lightblue","red","yellow","blue","green","grey","orange"))+
#   ggpubr::stat_compare_means(method="t.test",label="p.signif")+
#   scale_y_continuous(limits=c(0,0.4),expand = expansion(mult = c(0.05, 0.1)))+
#   theme_classic()+
#   ylab("Relative Abundance")+ggtitle("Baseline")+xlab("Cell Type")
# dev.off()

################ (2) After CIBERSORT comparison ################

cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/WANTED2_CIBERSORTx_Adjusted_TRM.txt",row.names = 1)
cibersort_res = cibersort_res[cell_annot$Patient,]
cibersort_res$response = ifelse(cell_annot$Response%in%c("Complete Response","CR","R","Partial Response","PR"),1,
                                ifelse(cell_annot$Response%in%c("Progressive Disease","PD","NR","SD"),0,"NE"))
cibersort_res$timepoint = cell_annot$timepoint
cibersort_res_long = melt(data.table::setDT(cibersort_res[,c(1,2,4,5,6,10,11)]),id.vars = c("response","timepoint"))
cibersort_res_long$response = as.factor(cibersort_res_long$response)
cibersort_res_long$timepoint = factor(cibersort_res_long$timepoint,levels=c("pre","post"))
cibersort_res_long = cibersort_res_long[which(cibersort_res_long$variable%in%c("TRM.T.cell","Exhausted.T.cell","TEM.T.cell","Naive.T.cell","Effector.T.cell")&cibersort_res_long$response!="NE"),]
cibersort_res_long$variable = factor(cibersort_res_long$variable,c("Exhausted.T.cell","TRM.T.cell","Effector.T.cell","TEM.T.cell","Naive.T.cell"))
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Figure5.NaiveEffec_RvsNR.pdf",
#     height=2.5,width=5)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.after_celltype_histogram.pdf",height=4,width=8)
ggplot(cibersort_res_long[which(cibersort_res_long$timepoint=="post"),],aes(x=response,y=value,fill=variable))+
  geom_boxplot(outlier.shape = NA)+
  # facet_wrap(~variable,ncol=5)+
  scale_fill_manual(values=c("blueviolet","orchid2","palegreen3","palegreen","lightblue"))+
  #ggpubr::stat_compare_means(method="t.test",label='p.signif')+
  scale_y_continuous(limits=c(0,0.4),expand = expansion(mult = c(0.05, 0.1)))+theme_classic()+
  #geom_jitter(position=position_dodge(),alpha=0.2)+
  ylab("Relative Abundance")+ggtitle("On Treatment")+xlab("Cell Type")
dev.off()

################ (3) Baseline some genes comparison ################

tmp = as.data.frame(
  cbind(c(t(allcount["GZMA",which(cell_annot$Response_binary==1&cell_annot$timepoint=="pre")]),
          t(allcount["GZMA",which(cell_annot$Response_binary==0&cell_annot$timepoint=="pre")])),
        c(t(allcount["FOXP3",which(cell_annot$Response_binary==1&cell_annot$timepoint=="pre")]),
          t(allcount["FOXP3",which(cell_annot$Response_binary==0&cell_annot$timepoint=="pre")])),
        c(rep("Responder",length(which(cell_annot$Response_binary==1&cell_annot$timepoint=="pre"))),
          rep("Non-responder",length(which(cell_annot$Response_binary==0&cell_annot$timepoint=="pre"))))))
tmp$V1 = as.numeric(tmp$V1)
tmp$V2 = as.numeric(tmp$V2)
names(tmp) = c("GZMA","FOXP3","Response")
tmp_long = melt(setDT(tmp),id.vars="Response")
ggplot(tmp_long,aes(x=Response,y=value,fill=Response))+geom_violin()+facet_wrap(~variable)+ggpubr::stat_compare_means()

################ (4) Baseline NK cells comparison ################

cibersort_res = read.delim("/Volumes/she4/hallmark/data/validation-data/WANTED1_CIBERSORTx_Adjusted_celltype.txt",row.names = 1)
cibersort_res = cibersort_res[cell_annot$Patient,]
cibersort_res$response = ifelse(cell_annot$Response%in%c("Complete Response","CR","Partial Response","PR","R"),1,
                                ifelse(cell_annot$Response%in%c("Progressive Disease","PD","NR","SD"),0,"NE"))
cibersort_res$timepoint = cell_annot$timepoint
cibersort_res_long = melt(data.table::setDT(cibersort_res[,c(1:22,26,27)]),id.vars = c("response","timepoint"))
cibersort_res_long$response = as.factor(cibersort_res_long$response)
cibersort_res_long$timepoint = factor(cibersort_res_long$timepoint,levels=c("pre","post"))
cibersort_res_long_select = cibersort_res_long[which(cibersort_res_long$variable%in%c("NK.cells.resting","NK.cells.activated")),]
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/Fig5.baseline_NK.pdf",height=3,width=6)
ggplot(cibersort_res_long_select[which(cibersort_res_long_select$response!="NE"&cibersort_res_long_select$timepoint=='pre'),],
       aes(x=response,y=value,fill=response))+
  geom_violin()+
  ggpubr::stat_compare_means(method="t.test",label="p.signif")+facet_wrap(~variable,scale="free_y")+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  ylab("Relative Abundance")+ggtitle("Baseline")+xlab("Response")+theme_classic()
dev.off()

# # after only:
# cell_annot_pre = cell_annot[which(cell_annot$timepoint%in%c("post")),]
# full_data = as.data.frame(cbind(ssGSEA[which(cell_annot$timepoint%in%c("post")),],cell_annot_pre))
# full_data = full_data[,c(names(full_data)[1:111],"Response_binary")]
# full_data = full_data[complete.cases(full_data),]
# full_data_long = melt(setDT(full_data),id.vars="Response_binary")
# full_data_long_MP = full_data_long[which(full_data_long$variable%in%names(gs_list_MP)),]
# full_data_long_MP$Response_binary = as.factor(full_data_long_MP$Response_binary)
# ggplot(full_data_long_MP,aes(x=Response_binary,y=value,fill=Response_binary))+geom_violin()+
#   facet_wrap(~variable)+ggpubr::stat_compare_means(method='t.test')+
#   stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",colour = "red")+
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
#


