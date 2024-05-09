library(msigdbr) #ok
library(org.Hs.eg.db) #ok
library(dplyr)
library(tidyr)
library(purrr)
#install.packages("readr")
library(readr)
library(stringr)
library(glue)
library(xml2)
#BiocManager::install("GEOquery")
library(GEOquery)
library("AnnotationDbi")
#BiocManager::install("biomaRt")
#library("hgu95av2.db")  
#library("biomaRt") 
#library(DESeq2)
#install.packages("NMF", force=T)
library(NMF)
library(ComplexHeatmap)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
#library(rowr)
library(xlsx)
library(qdapRegex)


path_data = "/Volumes/she4/hallmark/data/"
#path_data = "Z:/she4/hallmark/data/"

path_result = "/Volumes/she4/hallmark/results/"
#path_result = "Z:/she4/hallmark/results/"

# gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "C7"))
# result_MSig = matrix(0,length(unique(gene_sets$gs_name)),2)
# result_MSig[,1] = unique(gene_sets$gs_name)
# result_MSig = as.data.frame(result_MSig)
# for (p in 1:length(unique(gene_sets$gs_name))){
#   result_MSig$V2[p] = list(as.character(unique(gene_sets[which(gene_sets$gs_name==unique(gene_sets$gs_name)[p]),"entrez_gene"])))
# }
# result_MSig$origin = "MSigDB-C7"


# ------------------- Step 1: Download all C7 datasets

# Import the MSigDB XML file
mdb_doc <- read_xml(paste0(path_data,"msigdb_v7.5.1.xml.txt"))

# Extract the XML attributes and convert into a tibble
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/mdb_XML_description
# GENESET record attributes:
# * STANDARD_NAME: gene set name
# * SYSTEMATIC_NAME: gene set name for internal indexing purposes
# * CATEGORY_CODE: gene set collection code, e.g., C2
# * SUB_CATEGORY_CODE: gene set subcategory code, e.g., CGP
# * PMID: PubMed ID for the source publication
# * GEOID: GEO or ArrayExpress ID for the raw microarray data in GEO or ArrayExpress repository
# * EXACT_SOURCE: exact source of the set, usually a specific figure or table in the publication
# * GENESET_LISTING_URL: URL of the original source that listed the gene set members (all blank)
# * EXTERNAL_DETAILS_URL: URL of the original source page of the gene set
# * DESCRIPTION_BRIEF: brief description of the gene set
# * MEMBERS: list of gene set members as they originally appeared in the source
# * MEMBERS_SYMBOLIZED: list of gene set members in the form of human gene symbols
# * MEMBERS_EZID: list of gene set members in the form of human Entrez Gene IDs
# * MEMBERS_MAPPING: pipe-separated list of in the form of: MEMBERS, MEMBERS_SYMBOLIZED, MEMBERS_EZID


mdb_gs_ns <- xml_find_all(mdb_doc, xpath = ".//GENESET")
mdb_tbl <-
  tibble(
    gs_name = xml_attr(mdb_gs_ns, attr = "STANDARD_NAME"),
    gs_id = xml_attr(mdb_gs_ns, attr = "SYSTEMATIC_NAME"),
    gs_cat = xml_attr(mdb_gs_ns, attr = "CATEGORY_CODE"),
    gs_subcat = xml_attr(mdb_gs_ns, attr = "SUB_CATEGORY_CODE"),
    gs_pmid = xml_attr(mdb_gs_ns, attr = "PMID"),
    gs_geoid = xml_attr(mdb_gs_ns, attr = "GEOID"),
    gs_exact_source = xml_attr(mdb_gs_ns, attr = "EXACT_SOURCE"),
    gs_url = xml_attr(mdb_gs_ns, attr = "EXTERNAL_DETAILS_URL"),
    gs_description = xml_attr(mdb_gs_ns, attr = "DESCRIPTION_BRIEF"),
    gs_members = xml_attr(mdb_gs_ns, attr = "MEMBERS_MAPPING")
  ) %>%
  filter(gs_cat != "ARCHIVED")
mdb_tbl = mdb_tbl[which(mdb_tbl$gs_cat=="C7"),]
mdb_tbl = mdb_tbl[-which(mdb_tbl$gs_geoid==""),]
write.csv(mdb_tbl,paste0(path_data,"C7_source.csv"))

for (i in 1:length(unique(mdb_tbl$gs_geoid))){
  gse <- getGEO(unique(mdb_tbl$gs_geoid)[i], GSEMatrix = TRUE)
  #gse <- getGEOSuppFiles(unique(mdb_tbl$gs_geoid)[i])
  list = gse[[names(gse)[1]]]
  expression = list@assayData$exprs
  phenotype = list@phenoData@data
  feature = list@featureData@data
  rownames(expression)=feature$`Gene Symbol`
  
  source = gsub(" ", "", unique(phenotype$source_name_ch1))[1]
  if(grepl("/",source)){source = gsub("/","-",source)}
  
  if(unique(phenotype$organism_ch1)=="Homo sapiens"){
    write.csv(expression,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_expression_human.csv"))
    write.csv(phenotype,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_phenotype_human.csv"))
    write.csv(feature,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_feature_annotation_human.csv"))
  }else{
    write.csv(expression,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_expression_nonhuman.csv"))
    write.csv(phenotype,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_phenotype_nonhuman.csv"))
    write.csv(feature,paste0(path_data,"c7_GEO/",unique(mdb_tbl$gs_geoid)[i],"_",source,"_feature_annotation_nonhuman.csv"))
  }
}

# ------------------- Step 2: NMF for each dataset

# C7 gene sets
# gene_sets = as.data.frame(msigdbr::msigdbr(species = "Homo sapiens", category = "C7"))
# gene_sets = gene_sets[-which(gene_sets$gs_geoid==""),]

# NMF_program_all = list()

file_names <- dir(paste0(path_data,"c7_GEO"),full.names = TRUE,pattern = "Myeloid_expression_human.csv") #where you have your files

i=0
for (i in 1:length(file_names)){ # total 391 datasets: 47 human Myeloid
  
  geoid = qdapRegex::ex_between(file_names[i], "c7_GEO/", "_Myeloid")[[1]]
  
  expr = read.csv(paste0(path_data,"c7_GEO/",geoid,"_Myeloid_expression_human.csv"))
  pheno = read.csv(paste0(path_data,"c7_GEO/",geoid,"_Myeloid_phenotype_human.csv"))
  feature_annotation = read.csv(paste0(path_data,"c7_GEO/",geoid,"_Myeloid_feature_annotation_human.csv"))
  
  #expr$gene_symbol = feature_annotation[which(feature_annotation$X==expr$X),"Gene.Symbol"]
  expr$X = toupper(expr$X)
  expr = expr[!duplicated(expr$X),]
  
  if(length(which(is.na(expr$X=="")))==0){
    expr = expr
  }else{
    expr = expr[-which(expr$X==""),]
  }
  
  if(sum(is.na(expr$X))==0){
    rownames(expr) = toupper(expr$X)
  }else{
    expr = expr[-which(is.na(expr$X)),]
    rownames(expr) = toupper(expr$X)
  }
  expr = expr[,-grep("X",colnames(expr))]
  expr = expr[complete.cases(expr),]
  
  # Get rid of genes with 25% variance across samples
  expr = expr[-which(apply(expr,1,sd)<summary(apply(expr,1,sd))[2]),]
  # Get rid of RPS genes
  expr = expr[-grep("RPS",rownames(expr)),]
  expr = expr[-grep("RPL",rownames(expr)),]
  # Get rid of genes with weird names
  expr = expr[-grep("///",rownames(expr)),]
  
  numbers_only <- function(x) !grepl("\\D", x)

  if(numbers_only(rownames(expr)[1])){next}
  
  if(any(as.matrix(expr)<0)){next}
  
  # Extract the wanted genes from all gene sets that came out of this particular study
  #gs = gene_sets[which(gene_sets$gs_geoid==unique(mdb_tbl$gs_geoid)[i]),"gene_symbol"]
  
  # NMF
  #subset_expr = expr[which(rownames(expr)%in%gs),]
  subset_expr = expr

  # range = if(ncol(subset_expr)<10){
  #   seq(2,ncol(subset_expr))
  # }else{
  #   seq(2,10)  
  # } # need to find the optimal rank
  # rand = nmfEstimateRank(as.matrix(subset_expr),range=range, method='brunet', nrun=2, seed=123456)
  # #plot(rand)
  
  # take the top 50 genes from each component
  if(ncol(subset_expr)<4){next}
  i=i+1
  if(ncol(subset_expr)==4){
    decomp <- nmf(as.matrix(subset_expr),rank=4) 
    
    W = as.data.frame(decomp@fit@W)
    NMF_program = matrix(0,50,4)
    NMF_program_loadings = matrix(0,50,4) 
    for (r in 1:4){
      NMF_program[,r] = rownames(head(W[order(W[,r],decreasing = T),],50)) # top 50 genes from each PC
      NMF_program_loadings[,r] = head(W[order(W[,r],decreasing = T),],50)[,r]
    }
  }else if(ncol(subset_expr)==5){
    decomp <- nmf(as.matrix(subset_expr),rank=4) 
    decomp1 <- nmf(as.matrix(subset_expr),rank=5) 
    
    W = as.data.frame(decomp@fit@W)
    NMF_program = matrix(0,50,4)
    NMF_program_loadings = matrix(0,50,4) 
    for (r in 1:4){
      NMF_program[,r] = rownames(head(W[order(W[,r],decreasing = T),],50))
      NMF_program_loadings[,r] = head(W[order(W[,r],decreasing = T),],50)[,r]
    }

    W1 = as.data.frame(decomp1@fit@W)
    NMF_program1 = matrix(0,50,5)
    NMF_program_loadings1 = matrix(0,50,5) 
    for (r in 1:5){
      NMF_program1[,r] = rownames(head(W1[order(W1[,r],decreasing = T),],50))
      NMF_program_loadings1[,r] = head(W1[order(W1[,r],decreasing = T),],50)[,r]
    }
    NMF_program = cbind(NMF_program,NMF_program1)
    NMF_program_loadings = cbind(NMF_program_loadings,NMF_program_loadings1)
  }else{
    decomp <- nmf(as.matrix(subset_expr),rank=4) 
    decomp1 <- nmf(as.matrix(subset_expr),rank=5) 
    decomp2 <- nmf(as.matrix(subset_expr),rank=6) 
    
    W = as.data.frame(decomp@fit@W)
    NMF_program = matrix(0,50,4)
    NMF_program_loadings = matrix(0,50,4) 
    for (r in 1:4){
      NMF_program[,r] = rownames(head(W[order(W[,r],decreasing = T),],50))
      NMF_program_loadings[,r] = head(W[order(W[,r],decreasing = T),],50)[,r]
    }
    
    W1 = as.data.frame(decomp1@fit@W)
    NMF_program1 = matrix(0,50,5)
    NMF_program_loadings1 = matrix(0,50,5) 
    for (r in 1:5){
      NMF_program1[,r] = rownames(head(W1[order(W1[,r],decreasing = T),],50))
      NMF_program_loadings1[,r] = head(W1[order(W1[,r],decreasing = T),],50)[,r]
    }
    
    W2 = as.data.frame(decomp2@fit@W)
    NMF_program2 = matrix(0,50,6)
    NMF_program_loadings2 = matrix(0,50,6) 
    for (r in 1:6){
      NMF_program2[,r] = rownames(head(W2[order(W2[,r],decreasing = T),],50))
      NMF_program_loadings2[,r] = head(W2[order(W2[,r],decreasing = T),],50)[,r]
    }
    
    NMF_program = cbind(NMF_program,NMF_program1,NMF_program2)
    NMF_program_loadings = cbind(NMF_program_loadings,NMF_program_loadings1,NMF_program_loadings2)
  }
  #Heatmap(decomp@fit@H) # k by n

  # W = as.data.frame(decomp@fit@W)
  # NMF_program = matrix(0,50,range[which.max(rand$measures$cophenetic)])
  # for (r in 1:range[which.max(rand$measures$cophenetic)]){
  #   NMF_program[,r] = rownames(head(W[order(W[,r],decreasing = T),],50))
  # }
  
  write.csv(NMF_program,paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/humanMyeloid_genes_",unique(mdb_tbl$gs_geoid)[i],"_gene.csv"))
  write.csv(NMF_program_loadings,paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/humanMyeloid_genes_",unique(mdb_tbl$gs_geoid)[i],"_loadings.csv"))
  
  # NMF_program_all[[i]] = NMF_program
}
# problems: 47

# ------------------- Step 2.5: Define Robust NMFs
file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "gene.csv") #where you have your files
file_names_loadings <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_loadings.csv") #where you have your files
### Normalize the loading
for (i in 1:length(file_names_loadings)){
  sload = scale(read.csv(file_names_loadings[i],row.names = 1),scale = T,center = F)
  write.csv(sload,paste0(strsplit(file_names_loadings[i],"\\.")[[1]][1],"_scaled.csv"))
}
###
file_names_loadings <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_loadings_scaled.csv") #where you have your files
NMF_program_all <- do.call(cbind,lapply(file_names,read.csv))
NMF_program_all = NMF_program_all[,-grep("X",colnames(NMF_program_all))]
NMF_program_all_loadings <- do.call(cbind,lapply(file_names_loadings,read.csv))
NMF_program_all_loadings = NMF_program_all_loadings[,-grep("X",colnames(NMF_program_all_loadings))]

# (1) Robust across datasets: Keep NMF programs that had at least 20% similarity with any NMF program in any of the other datasets analyzed

for (i in 1:length(file_names)){
  df = read.csv(file_names[i])[,-1]
  df_loadings = read.csv(file_names_loadings[i])[,-1]
    
  index_rm = NULL # remove NMFs from this dataset
  for (c in 1:ncol(df)){
    index_rm = c(index_rm,unname(which(apply(NMF_program_all,2,function(x){sum(x==df[,c])})==50))) # remove itself, so that it's not compared against itself
  }
  NMF_program_all_rest = NMF_program_all[,-index_rm]
  NMF_program_all_loadings_rest = NMF_program_all_loadings[,-index_rm]
  
  rm = NULL
  for (c in 1:ncol(df)){
    if(sum(apply(NMF_program_all_rest,2,function(x){length(x=intersect(df[,c],x))/50})>0.2)==0){rm = c(rm,c)} # if no other NMFs overlap >20% with the family NMF, then remove
  }
  if(length(rm)==0){
    write.csv(df,paste0(strsplit(file_names[i],".csv")[1],"_robust.csv"))
    write.csv(df_loadings,paste0(strsplit(file_names_loadings[i],".csv")[1],"_robust.csv"))
  }else{
    write.csv(df[,-rm],paste0(strsplit(file_names[i],".csv")[1],"_robust.csv"))
    write.csv(df_loadings[,-rm],paste0(strsplit(file_names_loadings[i],".csv")[1],"_robust.csv"))
  }
}

# (2) Within each dataset, NMF programs were ranked by their Max similarity with NMFs from other datasets and 
##    selected in decreasing order. 
##    Once an NMF was selected, any other NMF within the dataset that had >= 20% overlap with the selected NMF was removed. 

file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_gene_robust.csv")
file_names_loadings <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_loadings_scaled_robust.csv")
index = NULL
for (i in 1:length(file_names)){
  if(sum(is.na(read.table(file_names[i],sep = ",")[,-1]))==50){index = c(index,i)}
  #if(is.null(ncol(read.table(file_names[i],sep = ",")[,-1]))){index = c(index,i)} # false
}
file_names = file_names[-index]
file_names_loadings = file_names_loadings[-index] 

NMF_program_all_robust <- do.call(cbind,lapply(file_names,read.csv))
NMF_program_all_robust = NMF_program_all_robust[,-grep("X",colnames(NMF_program_all_robust))]
NMF_program_all_loadings_robust <- do.call(cbind,lapply(file_names_loadings,read.csv))
NMF_program_all_loadings_robust = NMF_program_all_loadings_robust[,-grep("X",colnames(NMF_program_all_loadings_robust))]

for (i in 1:length(file_names)){
  
  df = read.csv(file_names[i])[,-1]
  df_loadings = read.csv(file_names_loadings[i])[,-1]
  
  # Scenario 1: If there's only one NMF in this dataset, then output this NMF directly
  if(is.null(ncol(df))){ 
    write.csv(df,paste0(strsplit(file_names[i],".csv")[1],"2.csv"))
    write.csv(df_loadings,paste0(strsplit(file_names_loadings[i],".csv")[1],"2.csv"));next
  }
  
  # Scenario 2: If only two NMFs left, compare the two, remove one or keep both
  if(ncol(df)==2){
    if(length(intersect(df[,1],df[,2]))/50>=0.2){
      redundancy_rm = rbinom(1,1,0.5)+1 #  remove either one
      df = df[,-redundancy_rm]
      df_loadings = df_loadings[,-redundancy_rm]
      write.csv(df,paste0(strsplit(file_names[i],".csv")[1],"2.csv"))
      write.csv(df_loadings,paste0(strsplit(file_names_loadings[i],".csv")[1],"2.csv"));next
    }else{
      df = df
      df_loadings = df_loadings
      write.csv(df,paste0(strsplit(file_names[i],".csv")[1],"2.csv"))
      write.csv(df_loadings,paste0(strsplit(file_names_loadings[i],".csv")[1],"2.csv"));next
    }
  }
  
  # Scenario 3: If more than two NMFs in the original data
  
  ## First, Remove the NMFs from this dataset itself 
  index_rm = NULL
  for (c in 1:ncol(df)){
    index_rm = c(index_rm,unname(which(apply(NMF_program_all_robust,2,function(x){sum(x==df[,c])})==50)))
  }
  NMF_program_all_robust_rest = NMF_program_all_robust[,-index_rm]
  NMF_program_all_robust_loadings_rest = NMF_program_all_loadings_robust[,-index_rm]
    
  ## Second, Iteratively choose the top 1 NMF, remove the NMFs with >=20% overlap in this dataset 
  keep = NULL
  keep_loadings = NULL
    while(!is.null(ncol(df))){ # do this until we have one column left or no column left # !is.null(ncol(df))|
      max = NULL
      for (c in 1:ncol(df)){ # this will give you a vector max with length ncol(df)
        max = c(max,max(apply(NMF_program_all_robust_rest,2,function(x){length(x=intersect(df[,c],x))/50})))
      }
      if(ncol(df)==2){ # if only two NMFs left 
        if(length(intersect(df[,-which.max(max)],df[,which.max(max)]))/50>=0.2){
          redundancy_rm = rbinom(1,1,0.5)+1
          df = df[,-redundancy_rm]
          df_loadings = df_loadings[,-redundancy_rm];break
        }else{
          df = df
          df_loadings = df_loadings;break
        }
      }
      if(length(which(apply(df[,-which.max(max)],2,function(x){length(intersect(df[,which.max(max)],x))/50})>=0.2))==0){ # nothing >= 20%, then keep iterate to the next top NMF
        keep = cbind(keep,df[,which.max(max)])
        keep_loadings = cbind(keep_loadings,df_loadings[,which.max(max)])
        df = df[,-which.max(max)]
        df_loadings = df_loadings[,-which.max(max)]
      }else{
        keep = cbind(keep,df[,which.max(max)])
        keep_loadings = cbind(keep_loadings,df_loadings[,which.max(max)])
        redundancy_rm = unname(which(apply(df[,-which.max(max)],2,function(x){length(x=intersect(df[,which.max(max)],x))/50})>=0.2)) # remove all NMFs with >= 20%
        df = df[,-c(which.max(max),redundancy_rm)] 
        df_loadings = df_loadings[,-c(which.max(max),redundancy_rm)]
      }

      if(is.null(ncol(df))|ncol(as.data.frame(df))==0){break}
    }
  
  write.csv(cbind(df,keep),paste0(strsplit(file_names[i],".csv")[1],"2.csv"))
  write.csv(cbind(df_loadings,keep_loadings),paste0(strsplit(file_names_loadings[i],".csv")[1],"2.csv"))
}


# ------------------- Step 3: Merge programs

file_names <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_gene_robust2.csv")
NMF_program_all <- do.call(cbind,lapply(file_names,read.csv))
NMF_program_all = NMF_program_all[,-grep("X",colnames(NMF_program_all))]

file_names_loadings <- dir(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_c7_NMF_programs/"),full.names = TRUE,pattern = "_loadings_scaled_robust2") # where you have your files
NMF_program_all_loadings <- do.call(cbind,lapply(file_names_loadings,read.csv))
NMF_program_all_loadings = NMF_program_all_loadings[,-grep("X",colnames(NMF_program_all_loadings))]


# jaccard_matrix = pbapply::pbsapply(1:ncol(NMF_program_all), FUN = Vectorize(function(a, b)(
#   length(intersect(unlist(NMF_program_all[,a]), unlist(NMF_program_all[,b])))/(length(unlist(NMF_program_all[,a])) + length(unlist(NMF_program_all[,b])) - length(intersect(unlist(NMF_program_all[,a]), unlist(NMF_program_all[,b]))))
# )))
# 
# overlap_matrix = outer(1:ncol(NMF_program_all), 1:ncol(NMF_program_all), FUN = Vectorize(function(a, b)(
#   length(intersect(unlist(NMF_program_all[,a]), unlist(NMF_program_all[,b]))))))
# any20 = function(x) sum(x<20)
# less_than_20_overlap_index = which(apply(overlap_matrix,2,any20)==(dim(overlap_matrix)[1]-1))

# Method 1: Hierarchical Clustering

# Problem: very uneven cluster size
# d = dist(jaccard_matrix)
# hc = fastcluster::hclust(d, method="average")
# memb <- cutree(hc, k = 50)
# venn.diagram(
#   x = list(NMF_program_all[,which(memb==1)[1]], NMF_program_all[,which(memb==1)[2]],
#            NMF_program_all[,which(memb==1)[3]], NMF_program_all[,which(memb==1)[4]],
#            NMF_program_all[,which(memb==1)[5]]),
#   category.names = c("Set 1" , "Set 2 ","Set 3" , "Set 4 ", "Set 5 "),
#   filename = '~/desktop/venn.png',
#   output=TRUE
# )

# Method 2: Cluster base to be NMFs with top 50 Jaccard distance

# temp = apply(jaccard_matrix,2,sum)/2
# top = head(sort(temp,decreasing = T),50)
# cluster_base = which(temp%in%top)
# 
# NMF_program_all_minus_base = NMF_program_all[,-cluster_base]
# 
# result = matrix(0,length(cluster_base),ncol(NMF_program_all_minus_base))
# for (b in 1:length(cluster_base)){
#   for (i in 1:ncol(NMF_program_all_minus_base)){
#     result[b,i] = length(intersect(unlist(NMF_program_all_minus_base[,i]), 
#                                    unlist(NMF_program_all[,cluster_base[b]])))
#   }
# }
# tmp = apply(result,2,which.max)

# Method 3: Tirosh Method

to_delete = NULL
to_add = NULL

# if(length(to_delete) == 0){
#   NMF_program_all = NMF_program_all
#   NMF_program_all_loadings = NMF_program_all_loadings
# }else{
#   NMF_program_all = NMF_program_all[,-to_delete]
#   NMF_program_all_loadings = NMF_program_all_loadings[,-to_delete]
# }

iter = 0

while(!is.null(ncol(NMF_program_all))){
  
  iter = iter + 1 
  
  # calculate pairwise overlap
  overlap_matrix = outer(1:ncol(NMF_program_all), 1:ncol(NMF_program_all), FUN = Vectorize(function(a, b)(
    length(intersect(unlist(NMF_program_all[,a]), unlist(NMF_program_all[,b]))))))
  
  # overlap_matrix = pbapply::pbsapply(1:ncol(NMF_program_all), FUN = Vectorize(function(a, b)(
  #   length(intersect(unlist(NMF_program_all[,a]), unlist(NMF_program_all[,b]))))))
  
  sum_overlap = apply(overlap_matrix,2,sum)
  
  # top1 = NULL
  
  top1 = which.max(sum_overlap) # find the top 1 NMF 
  if(sum(overlap_matrix[,tail(top1)]>10)>5){ # if # of overlapping NMF (defined as >20% overlap) > 5
    #to_add = which(overlap_matrix[-tail(top1),tail(top1)]==sort(overlap_matrix[-tail(top1),tail(top1)],decreasing = T)[1])
    to_add = c(to_add,which(overlap_matrix[-tail(top1),tail(top1)]==sort(overlap_matrix[-tail(top1),tail(top1)],decreasing = T)[1]))
    to_delete = c(top1,to_add)
    # df = as.data.frame(cbind(NMF_program_all[,c(to_delete)],NMF_program_all_loadings[,c(to_delete)]))
    df = as.data.frame(cbind(unname(unlist(NMF_program_all[,c(to_delete)])),
                             unname(unlist(NMF_program_all_loadings[,c(to_delete)]))))
    #
    intersect = intersect(unname(unlist(NMF_program_all[,c(top1)])),unname(unlist(NMF_program_all[,c(to_add)]))) # find intersect genes between top1 and added one
    stretch_num = 50 - length(intersect)
    
    MP_data1 =  df[which(df$V1%in%intersect),] 
    df_nointersect = df[-which(df$V1%in%intersect),]
    stretch = head(df_nointersect[order(df_nointersect$V2,decreasing = T),1],stretch_num)
    MP_data2 =  df[which(df$V1%in%stretch),]
    MP = c(intersect,stretch)
    MP_data = rbind(MP_data1,MP_data2)
    MP_data = MP_data[!duplicated(MP_data$V1),]
    #
    # MP_data =  head(df[order(df$V2,decreasing = T),],50)
    # MP = head(df[order(df$V2,decreasing = T),1],50) # update MP for the first time: here we could cap the # of the genes in MP, take the union?
    
    overlap = NULL
    for (i in 1:ncol(NMF_program_all[,-to_delete])){
      overlap = c(overlap,length(intersect(MP,NMF_program_all[,-to_delete][,i])))
    }
    if(!any(overlap>10)){write.csv(MP_data,paste0(path_result,"c7_NMF_programs/filter_genes_cluster_methodTirosh2/","MP_",top1,"_iter_",iter,"_based.csv"))}
    
    while(any(overlap>10)){
      filter_genes_to_add = unname(which(apply(NMF_program_all, 2, identical, NMF_program_all[,-to_delete][,which.max(overlap)]))) # add another NMF
      # filter_genes_to_add = NMF_program_all[,-to_delete][,which.max(overlap)]
      # NMF_program_all[,-delete][,which.max(overlap)]
      to_delete = c(to_delete,filter_genes_to_add)
      df = as.data.frame(rbind(MP_data, # MP
                               cbind(unname(unlist(NMF_program_all[,filter_genes_to_add])),unname(unlist(NMF_program_all_loadings[,filter_genes_to_add])))))
      #
      intersect = intersect(MP_data$V1,unname(unlist(NMF_program_all[,c(filter_genes_to_add)])))
      stretch_num = 50 - length(intersect)
      df_nointersect = df[-which(df$V1%in%intersect),]
      stretch = head(df_nointersect[order(df_nointersect$V2,decreasing = T),1],stretch_num)
      MP = c(intersect,stretch)
      MP_data1 =  df[which(df$V1%in%intersect),]
      MP_data2 =  df[which(df$V1%in%stretch),]
      MP_data = rbind(MP_data1,MP_data2)
      MP_data = MP_data[!duplicated(MP_data$V1),]
      #
      # MP_data = head(df[order(df$V2,decreasing = T),],50)
      # MP = head(df[order(df$V2,decreasing = T),1],50) # update MP again
      
      overlap = NULL
      for (i in 1:ncol(NMF_program_all[,-to_delete])){
        overlap = c(overlap,length(intersect(MP,NMF_program_all[,-to_delete][,i])))
      }
    }
    write.csv(MP_data,paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/","scaled_MP_",top1,"_iter_",iter,"_based.csv"))
  }else{
    to_delete = c(to_delete,top1)
    write.csv(cbind(NMF_program_all[,top1],NMF_program_all_loadings[,top1]),paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/","scaled_MP_",top1,"_iter_",iter,"_based.csv"))
    # end clustering
  }
  NMF_program_all = NMF_program_all[,-to_delete]
  NMF_program_all_loadings = NMF_program_all_loadings[,-to_delete]
  to_delete = NULL
  to_add = NULL
}


# ------------------- Step 4: Annotate programs

# Following method 2
# for (b in 1:length(names(table(tmp)))){
#   cols <- c("SYMBOL", "GENENAME","ENTREZID")
#   i=names(table(tmp))[b]
#   map = AnnotationDbi::select(org.Hs.eg.db, keys=unname(unlist(as.vector(NMF_program_all[,c(cluster_base[as.numeric(i)],which(as.character(tmp)==i))]))), 
#                               columns=cols, keytype="SYMBOL")
#   map = map[complete.cases(map),]
#   map = map[-which(duplicated(map$SYMBOL)),]
#   KEGG = enrichKEGG(gene = map$ENTREZID)
#   edox2 <- pairwise_termsim(KEGG)
#   p1 <- treeplot(edox2)
#   pdf(paste0(path_result,"c7_NMF_programs/","cluster_",b,".pdf"),width=15,height = 8)
#   print(p1)
#   dev.off()
#   
#   #For b=18/24/37
#   # pdf(paste0(path_result,"c7_NMF_programs/","cluster_",37,".pdf"),width=8,height = 8)
#   # print(barplot(KEGG))
#   # dev.off()
# }
# edox <- setReadable(KEGG, 'org.Hs.eg.db', 'ENTREZID')
# cnetplot(edox, showCategory = 5, categorySize="pvalue", colorEdge = TRUE,node_label="gene")

# Following method 3

# msigdbr_df_c3 <- msigdbr(species = "human", category = "C3")
# msigdbr_df_c3 = as.data.frame(msigdbr_df_c3)
# msigdbr_df_c3 = msigdbr_df_c3[grep("TFT:",msigdbr_df_c3$gs_subcat),]
# term2gene_c3 <- msigdbr_df_c3[, c("gs_name", "entrez_gene")]
# term2name_c3 <- msigdbr_df_c3[, c("gs_name", "gs_description")]

file_names <- dir(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2"),full.names = TRUE,pattern = "scaled_MP")[1:19] #where you have your files

tmp_MF=list()
tmp_BP=list()
tmp_KEGG=list()
tmp_H=list()

tmp = list()
enrich_result = list()

install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method","auto")

msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
t2g = cbind(rep(names(pathwaysH)[1],length(pathwaysH[[1]])),pathwaysH[[1]])
for (i in 1:50){
  t2g = rbind(t2g,cbind(rep(names(pathwaysH)[i],length(pathwaysH[[i]])),pathwaysH[[i]]))
}


# length(file_names)
for (i in 1:length(file_names)){
  cols <- c("SYMBOL", "GENENAME","ENTREZID")
  #i=names(table(tmp))[b]
  map = AnnotationDbi::select(org.Hs.eg.db, 
                              keys=read.csv(file_names[i])[,2], 
                              columns=cols, keytype="SYMBOL")
  map = map[complete.cases(map),]
  # map = map[-which(duplicated(map$SYMBOL)),]
  GO = enrichGO(gene = map$ENTREZID,"org.Hs.eg.db",ont="MF")
  if(length(GO@result[which(GO@result$p.adjust<0.05),2]) == 0){
    tmp_MF[[i]] = c(NA,NA,NA)
  }else{
    GO@result = GO@result[order(GO@result$p.adjust,decreasing = F),]
    tmp_MF[[i]] = cbind(paste("MF:",GO@result[which(GO@result$p.adjust<0.05),2]),
                        as.numeric(GO@result[which(GO@result$p.adjust<0.05),6]),
                        as.numeric(GO@result[which(GO@result$p.adjust<0.05),9]))
  }

  GO2 = enrichGO(gene = map$ENTREZID,"org.Hs.eg.db",ont="BP")
  if(length(GO2@result[which(GO2@result$p.adjust<0.05),2]) == 0){
    tmp_BP[[i]] = c(NA,NA,NA)
  }else{
    GO2@result = GO2@result[order(GO2@result$p.adjust,decreasing = F),]
    tmp_BP[[i]] = cbind(paste("BP:",GO2@result[which(GO2@result$p.adjust<0.05),2]),
                        as.numeric(GO2@result[which(GO2@result$p.adjust<0.05),6]),
                        as.numeric(GO2@result[which(GO2@result$p.adjust<0.05),9]))
  }

  KEGG = enrichKEGG(gene = map$ENTREZID,organism = "hsa")
  if(length(KEGG@result[which(KEGG@result$p.adjust<0.05),2]) == 0){
    tmp_KEGG[[i]] = c(NA,NA,NA)
  }else{
    KEGG@result = KEGG@result[order(KEGG@result$p.adjust,decreasing = F),]
    tmp_KEGG[[i]] = cbind(paste("KEGG:",KEGG@result[which(KEGG@result$p.adjust<0.05),2]),
                          as.numeric(KEGG@result[which(KEGG@result$p.adjust<0.05),6]),
                          as.numeric(KEGG@result[which(KEGG@result$p.adjust<0.05),9]))
  }

  H = clusterProfiler::enricher(map$ENTREZID,TERM2GENE = t2g)
  if(length(H@result[which(H@result$p.adjust<0.05),2]) == 0){
    tmp_H[[i]] = c(NA,NA,NA)
  }else{
    H@result = H@result[order(H@result$p.adjust,decreasing = F),]
    tmp_H[[i]] = cbind(paste("H:",H@result[which(H@result$p.adjust<0.05),2]),
                       as.numeric(H@result[which(H@result$p.adjust<0.05),6]),
                       as.numeric(H@result[which(H@result$p.adjust<0.05),9]))
  }
  # C3 <- enricher(map$ENTREZID, TERM2GENE=term2gene_c3,TERM2NAME=term2name_c3)
  # tab.dis <- as.data.frame(C3)
  # if(length(subset(tab.dis, Count>5)[,2]) == 0){
  #   tmp_C3[[i]] = c(NA,NA,NA)
  # }else{
  #   tmp_C3[[i]] <- cbind(paste("C3-TFT:",subset(tab.dis, Count>5)[,2]),
  #                        as.numeric(subset(tab.dis, Count>5)[,6]),
  #                        as.numeric(subset(tab.dis, Count>5)[,9]))
  # }
  cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
      rbind(x, matrix(, n-nrow(x), ncol(x))))) 
  }
  # res = as.data.frame(cbind.fill(read.csv(file_names[i])[,2],rbind(tmp_BP[[i]],tmp_MF[[i]],tmp_KEGG[[i]]),fill = NA))
  # colnames(res) = c("Genes","Enriched Pathways","Adj_p","Counts",'')
  # res$Adj_p = as.numeric(res$Adj_p)
  # res$Counts = as.numeric(res$Counts)
  # write_xlsx(res, file=paste0(path_result,"c7_NMF_programs/filter_genes_cluster_methodTirosh2/res.xlsx"), sheetName=paste0("MP",i), append=TRUE, row.names=FALSE)

  
  tmp[[i]] = unlist(c(ifelse(is.null(dim(tmp_MF[[i]])),"",list(tmp_MF[[i]][,1])),
               ifelse(is.null(dim(tmp_BP[[i]])),"",list(tmp_BP[[i]][,1])),
               ifelse(is.null(dim(tmp_KEGG[[i]])),"",list(tmp_KEGG[[i]][,1])),
               ifelse(is.null(dim(tmp_H[[i]])),"",list(tmp_H[[i]][,1]))))
  saveRDS(tmp,"/Volumes/she4/hallmark/tmp.rds")
  
  df = rbind(GO@result[which(GO@result$p.adjust<0.05),],
             GO2@result[which(GO2@result$p.adjust<0.05),],
             KEGG@result[which(KEGG@result$p.adjust<0.05),],
             H@result[which(H@result$p.adjust<0.05),])
  if(nrow(df)==0){enrich_result[[i]] = NULL;next}
  for (r in 1:nrow(df)){
    df$geneID[r] = paste0(AnnotationDbi::select(org.Hs.eg.db, 
                                                keys=unlist(str_split(df$geneID[r],"/")), 
                                                columns=cols, keytype="ENTREZID")$SYMBOL,collapse = "/")
  }
  enrich_result[[i]] = df
  saveRDS(enrich_result,paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_enrich_result.rds"))
}

res = do.call(cbind.fill,tmp)
for(i in 1:ncol(res)){
  res[,i][duplicated(res[,i])] <- NA
}
write.csv(res,paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_res_new.csv"))

MPs = NULL
for (i in 1:length(file_names)){
  MPs = cbind(MPs,as.matrix(read.csv(file_names[i],row.names = 1)))
}
write.csv(MPs,paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_MPs_weight.csv"))

MPs = NULL
for (i in 1:length(file_names)){
  MPs = cbind(MPs,as.matrix(read.csv(file_names[i],row.names = 1)[,1]))
}
write.csv(MPs,paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_MPs.csv"))

enriched_M = readRDS(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/scaled_enrich_result.rds"))
enriched_L = readRDS(paste0(path_result,"c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/scaled_enrich_result.rds"))

library(xlsx)
for (i in 1:9){
  write.xlsx(enriched_M[[i]],
             "/Volumes/she4/hallmark/results/c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/enriched_M.xlsx",
             sheetName = as.character(i),append = T)
}
for (i in 1:19){
  write.xlsx(enriched_L[[i]],
             "/Volumes/she4/hallmark/results/c7_NMF_programs/humanLymphoid_genes_cluster_methodTirosh2/enriched_L.xlsx",
             sheetName = as.character(i),append = T)
}

jaccard_matrix = outer(1:9,1:9,FUN = Vectorize(function(a, b)(
  length(intersect(unlist(MPs[,a]), unlist(MPs[,b]))))))
mybreak=c(0,0.5,1)
mycolor = c("white","salmon","violet")
col_fun = colorRamp2(mybreak, mycolor)
rownames(jaccard_matrix) = paste0("MP",seq(1,9,1))
colnames(jaccard_matrix) = paste0("MP",seq(1,9,1))
pdf(paste0(path_result,"c7_NMF_programs/humanMyeloid_genes_cluster_methodTirosh2/MPs_overlap.pdf"))
Heatmap(jaccard_matrix/50,name = "Percentage overlap",col=col_fun,
        show_column_dend = F,show_row_dend = F,
        show_column_names = T,show_row_names = T,
        column_title = "Overlap among 9 Myeloid MPs")
dev.off()
