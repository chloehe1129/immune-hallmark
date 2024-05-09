library(rentrez)
library(GEOquery)
library(dplyr)
library(ggplot2)

mdb_tbl = read.csv(paste0(path_data,"C7_source.csv"),row.names = 1)

res = matrix(0,length(unique(mdb_tbl$gs_geoid)),1+1+1+1+1)
colnames(res) = c("journal","organism","immune_type","samples","genes")
rownames(res) = unique(mdb_tbl$gs_geoid) # 104
for (i in 105:length(unique(mdb_tbl$gs_geoid))){
  gse <- getGEO(unique(mdb_tbl$gs_geoid)[i], GSEMatrix = TRUE)
  #gse <- getGEOSuppFiles(unique(mdb_tbl$gs_geoid)[i])
  list = gse[[names(gse)[1]]]
  expression = list@assayData$exprs
  phenotype = list@phenoData@data
  feature = list@featureData@data

  pmid <- as.vector(unique(mdb_tbl[which(mdb_tbl$gs_geoid==unique(mdb_tbl$gs_geoid)[i]),"gs_pmid"]))[[1]]
  record <- entrez_fetch(db="pubmed", id=pmid, rettype="xml")

  res[i,1] = strsplit(strsplit(record,"<Title>")[[1]][2],"</Title>")[[1]][1]
  res[i,2] = unique(phenotype$organism_ch1)
  res[i,3] = gsub(" ", "", unique(phenotype$source_name_ch1))[1]
  res[i,4] = ncol(expression)
  res[i,5] = nrow(expression)
  print(i)
}
write.csv(res,"/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/res.csv")

res = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/res.csv",row.names = 1)
df = as.data.frame(cbind(c("Homo sapiens (35%)","Mus musculus (65%)"),round(c(sum(res[,2]=="Homo sapiens")/nrow(res),sum(res[,2]!="Homo sapiens")/nrow(res)),2)))
names(df) = c("Organism","Proportion")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/pie.pdf")
ggplot(df, aes(x="", y=Proportion, fill=Organism)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + xlab("")+ylab("")+theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.text = element_text(size=15))
dev.off()

res_homo = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/homo.csv",row.names = 1)
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/homo_samples.pdf")
hist(as.numeric(res_homo[,5]),main="Distribution of sample sizes across Human datasets",col="lightskyblue",xlab = "Sample size",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
text(x=80,y=30,labels = paste0("Mean: ",round(summary(res_homo[,5])[4],0)),cex=1.5)
text(x=80,y=25,labels = paste0("Quantile: [",round(summary(res_homo[,5])[2],0),",",round(summary(res_homo[,5])[5],0),"]"),cex=1.5)
dev.off()
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/homo_genes.pdf")
hist(as.numeric(res_homo[,6]),main="Distribution of feature count across Human datasets",col="lightskyblue",xlab = "Feature counts",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
text(x=500000,y=30,labels = paste0("Mean: ",round(summary(res_homo[,6])[4],0)),cex=1.5)
text(x=500000,y=25,labels = paste0("Quantile: [",round(summary(res_homo[,6])[2],0),",",round(summary(res_homo[,6])[5],0),"]"),cex=1.5)
dev.off()

# res_mus = res[which(res[,2]!="Homo sapiens"),]
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/mus_samples.pdf")
# hist(as.numeric(res_mus[,4]),main="Distribution of sample sizes across Mouse datasets",col="lavender",xlab = "Sample size")
# dev.off()
# pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/mus_genes.pdf")
# hist(as.numeric(res_mus[,5]),main="Distribution of feature count across Human datasets",col="lavender",xlab = "Feature counts")
# dev.off()

res1 = as.data.frame(res)

res1 <- res1 %>%
  group_by(journal) %>%
  summarise(count = n()) %>%
  mutate(per=count/sum(count)) %>% 
  ungroup()

df = as.data.frame(cbind(ifelse(as.data.frame(res1)$per>0.05,as.data.frame(res1)$journal,"Others"),res1$per))
names(df) = c("Journal","Proportion")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/pie_journal.pdf")
ggplot(df, aes(x="", y=Proportion, fill=Journal)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + xlab("")+ylab("")+theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.text = element_text(size=13))
dev.off()

res2 = read.csv("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/homo.csv")
res2 <- res2  %>%
  group_by(X.1) %>%
  summarise(count = n()) %>%
  mutate(per=count/sum(count)) %>% 
  ungroup()
labels = as.data.frame(res2)$X.1
labels[1] = "Others"

df = as.data.frame(cbind(labels,res2$per))
names(df) = c("CellType","Proportion")
pdf("/Users/she4/Library/CloudStorage/OneDrive-InsideMDAnderson/ken/hallmark/manuscript/supplement_figures/celltype_pie.pdf")
ggplot(df, aes(x="", y=Proportion, fill=CellType)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + xlab("")+ylab("")+theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        legend.text = element_text(size=13))
dev.off()

