options(stringsAsFactors=F)
library(DESeq2)
library(data.table)
library(dplyr)


Cond_fordeg<-read.table("/home/ha7477/share/to_imgkono2/result/siLEF1_ATAC/deg/meta.txt",header=T,sep="\t")
Count_fordeg<-fread("/home/ha7477/share/to_imgkono2/result/siLEF1_ATAC/featurecounts_summary/siLEF1_ATAC_all_featurecounts.txt",header=T,data.table=F,sep="\t",check.names = F)

rownames(Count_fordeg)<-Count_fordeg$peak
Count_fordeg<-Count_fordeg[,-1]
head(Count_fordeg,n=2)


colD<-Cond_fordeg[,c("sample_ID","condition","time","Donor")]
colnames(colD)<-c("sample_id","si","hr","donor")
colD$si<-factor(colD$si,levels=c("siCT","siLEF1"))
colD$hr<-factor(colD$hr,levels=c("24hr","48hr","72hr"))
colD$donor<-factor(colD$donor)

head(colD)

###permutation###
#tar_Count2<-tar_Count[,sample(1:ncol(tar_Count),ncol(tar_Count),replace = F)]
#colnames(tar_Count2)<-colnames(tar_Count)
#tar_Count<-tar_Count2
##################

dds <- DESeqDataSetFromMatrix(countData = Count_fordeg, colData = colD, design = ~ si + hr +donor )

#peak count filter
keep <- rowSums(counts(dds)) >= ncol(counts(dds))*10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds, contrast = c("si","siLEF1","siCT"))

res2<-res[which(res$padj<0.05),]
dim(res2)

