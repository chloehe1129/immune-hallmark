library(dplyr)
library(Seurat)

path_data = "/Volumes/she4/hallmark/data/"
path_result = "/Volumes/she4/hallmark/results/"

path_data = "/rsrch4/home/bcb/she4/hallmark/data/"
path_result = "/rsrch4/home/bcb/she4/hallmark/results/"

# Basic function to convert mouse to human gene names
convert_mouse_to_human <- function(gene_list) { 
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}

all_dic = list.files(paste0(path_data,"cytokine-dictionary/"),full.names = T)
all_dic_names = list.files(paste0(path_data,"cytokine-dictionary/"))

for(i in 1:length(all_dic)){
  dic = readRDS(all_dic[i])
  exp = dic@assays$RNA@counts
  gene_convert = convert_mouse_to_human(rownames(dic))
  gene_convert = gene_convert[!duplicated(gene_convert[,2]),]
  exp = exp[gene_convert[,1],]
  rownames(exp) = gene_convert[,2]
  meta = dic@meta.data
  dic = CreateSeuratObject(exp,meta.data=meta)
  saveRDS(dic,paste0(path_data,"cytokine-dictionary/human_",all_dic_names[i]))
}

all_human_dic = list.files(paste0(path_data,"cytokine-dictionary/"),full.names = T,pattern="human")
dic = readRDS(all_human_dic[13])
table(dic$sample)

a = readRDS(all_dic[31])
table(a$sample)
