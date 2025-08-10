rm(list=ls())
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(DESeq2)
library(org.Mm.eg.db) 
library(openxlsx)
library(dplyr)
library(readr)

#Reading in the gene expression matrix for samples from different groups.
#repeats.count.list like :DT_HFSC_a1	/Path/to/DT_HFSC_a1.repeats.count.simple.txt
setwd("/public/home/changying/2025_project/bulkRNA_2025/mouse_bulk_RNA_cuijun_2025_3_21/example/output/4.DEG/")
FeatureCountsList <- read.table("/public/home/changying/2025_project/bulkRNA_2025/mouse_bulk_RNA_cuijun_2025_3_21/example/code/featurecount.file.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(FeatureCountsList) <- c("sample","pathway") 
countslist <- list()  
for(i in 1:8){
  print(i)
  CountData <- lapply(FeatureCountsList$pathway[i], read.delim, header = TRUE, row.names = "Geneid")
  CountData <- do.call(cbind, CountData)
  colnames(CountData)[2] <- FeatureCountsList$sample[i]
  countslist[[i]] <- CountData
}

#Integrating ERV expression matrices from samples of different groups.
Matrix_erv<-merge(countslist[[1]],list(countslist[[2]],countslist[[3]],countslist[[4]]),by = "row.names", all = T)
rownames(Matrix_erv) <- Matrix_erv[,1]
Matrix_erv <- Matrix_erv[,-1]
Matrix_erv <- Matrix_erv[, !grepl("Length", names(Matrix_erv))]
colnames(Matrix_erv) <- c("6OHDA_DT_1_HFSC", "6OHDA_DT_3_HFSC","PBS_HFSC_a11","PBS_HFSC_a12" )
write.csv(Matrix_erv, file = "erv_row_Matrix_data.csv", sep = ",", row.names = TRUE)

#Integrating gene expression matrices from samples of different groups.
Matrix_mRNA<-merge(countslist[[5]],list(countslist[[6]],countslist[[7]],countslist[[8]]),by = "row.names", all = T)
rownames(Matrix_mRNA) <- Matrix_mRNA[,1]
Matrix_mRNA <- Matrix_mRNA[,-1]
Matrix_mRNA <- Matrix_mRNA[, !grepl("Length", names(Matrix_mRNA))]
colnames(Matrix_mRNA) <- c("6OHDA_DT_1_HFSC", "6OHDA_DT_3_HFSC","PBS_HFSC_a11","PBS_HFSC_a12" )
write.csv(Matrix_mRNA, file = "mRNA_row_Matrix_data.csv", sep = ",", row.names = TRUE)

#Differential ERV Expression Analysis between Different Groups by DESeq2
data1 <- Matrix_erv
condition1 <- factor(c(rep("group1",2),rep("group2",2)))
colData1 <- data.frame(row.names=colnames(data1), condition1) 
dds1 <- DESeqDataSetFromMatrix(countData = data1,colData = colData1,design = ~condition1)
dds1 <- DESeq(dds1)
res1 <- list()
res1[[1]] <- results(dds1,contrast=c("condition1","group2","group1")) 
#res1[[2]] <- results(dds,contrast=c("condition1","group2","group3"))

#load annotation file Mmus38.txt
gtf_in <- read_delim("/public/home/changying/project/mouse_bulk_RNA_cuijun_2024-2-1/Cell_2021_Lima-Junior/ref/Mmus38.txt",
                     "\t", escape_double = FALSE, col_names = FALSE,
                     comment = "#", trim_ws = TRUE)

colnames(gtf_in) <- c("ID","chr","start","end","strand","AA_length","method","N_cnt","MetID","Met_AA_length","HMM_profile","Viral_BLAST","NR_BLAST","EVE_BLAST","RetroTector","Repbase")
gtf_in <- gtf_in[-1,]
#Differential ERV annotation and saving differential comparison result files.
filename <- c("group2vsgroup1")
for (i in 1){
    resOrdered <- res1[[i]][order(-res1[[i]]$log2FoldChange),]
    resOrdered <- as.data.frame(resOrdered)
    annotation1 <- gtf_in[match(rownames(resOrdered), gtf_in$ID), ]
    res_ann <- cbind(resOrdered, annotation1[1:16])
    file_pre <- filename[i]
    file <-paste0(file_pre, "_DEseq2_erv.csv")
    write.csv(res_ann, file = file, sep = ",", row.names = TRUE)
}

#Differential gene Expression Analysis between Different Groups by DESeq2
data2 <- Matrix_mRNA
dds2 <- DESeqDataSetFromMatrix(countData = data2,colData = colData1,design = ~condition1)
dds2 <- DESeq(dds2)
res <- list()
res[[1]] <- results(dds2,contrast=c("condition1","group2","group1")) 
#res[[2]] <- results(dds2,contrast=c("condition1","group2","group3"))

#Differential gene annotation and saving differential comparison result files.
for (i in 1){
    res[[i]]$change = ifelse(res[[i]]$pvalue < 0.05 & abs(res[[i]]$log2FoldChange) >= 1, ifelse(res[[i]]$log2FoldChange > 1,'Up','Down'), 'Stable')
    resOrdered <- res[[i]][order(res[[i]]$padj),]
    resOrdered <- as.data.frame(resOrdered)
    file_pre <- filename[i]
    file <-paste0(file_pre, "_DEseq2_mRNA.csv")
    
    write.csv(resOrdered, file = file, sep = ",", row.names = TRUE)
}

