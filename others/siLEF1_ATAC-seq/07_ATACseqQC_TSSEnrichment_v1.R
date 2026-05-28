args   <- commandArgs(trailingOnly = T) 

bamfile=args[1]

options(stringsAsFactors=F)
library(ATACseqQC)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)


bamfile.labels <- gsub(".bam", "", basename(bamfile))


possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))

bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags

gal <- readBamFile(bamfile, tag=tags, asMates=TRUE, bigFile=FALSE)
list("n")[[1]]
gal1 <- shiftGAlignmentsList(gal)

tsse <- TSSEscore(gal1, txs)
TSSenrichmentScore<-tsse$TSSEscore

write.table(data.frame(t(c(bamfile.labels,TSSenrichmentScore))),paste0(bamfile,".TSSenrichmentscore"),quote=F,append=F,row.names=F,col.names=F)

pdf(paste0(bamfile,"_fragsizedist.pdf"))
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()

