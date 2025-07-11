options(stringsAsFactors=F)
library(DESeq2)
library(data.table)
library(dplyr)


Cond_fordeg<-read.table("/home/ha7477/share/to_imgkono2/result/siLEF1_RNA/deg/meta.txt",header=T,sep="\t")
Cond_fordeg$Sample<-paste0("FC01127_",Cond_fordeg$sample_ID)

su2<-fread("/home/ha7477/share/to_imgkono2/result/siLEF1_RNA/summary/count_rnaseqc2_summary.txt",data.table=F,header=T,check.names = F)
Count_fordeg<-su2[,-c(1,2)]
rownames(Count_fordeg)<-su2[,1]




colD<-Cond_fordeg[,c("Sample","condition","time","Donor")]
colnames(colD)<-c("sample_id","si","hr","donor")
colD$si<-factor(colD$si,levels=c("siCT","siLEF1"))
colD$hr<-factor(colD$hr,levels=c("24hr","48hr","72hr"))
colD$donor<-factor(colD$donor)

colD

###permutation###
#tar_Count2<-tar_Count[,sample(1:ncol(tar_Count),ncol(tar_Count),replace = F)]
#colnames(tar_Count2)<-colnames(tar_Count)
#tar_Count<-tar_Count2
##################

dds <- DESeqDataSetFromMatrix(countData = Count_fordeg, colData = colD, design = ~ si + hr +donor )

keep <- rowSums(counts(dds)) >= ncol(counts(dds))*10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds, contrast = c("si","siLEF1","siCT"))

res2<-res[which(res$padj<0.05),]
dim(res)
dim(res2)

