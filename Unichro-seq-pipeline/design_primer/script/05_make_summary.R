args <- commandArgs(trailingOnly = T)
IFILE <- as.character(args[1])
ODIR <- as.character(args[2])
   # IFILE="target/demo/demo.input"; ODIR="target/demo"

 #original input file
input <- read.table(IFILE, header=F)
colnames(input) <- c("Name", "CRE", "Strand_Fw1_Rv2","Dist_snp")

 #primer info
primer <- read.table(paste0(ODIR,"/target.primers"), header=T, sep="\t")

 #map info using all primer seq
map <- read.table(paste0(ODIR,"/g38_map.allseq.summary.txt"), header=T, sep="\t")
map$map_info <- paste0(map$Chr, ":", map$Position, ":", map$Strand)
allid <- as.character(unique(map$ID))
res <- data.frame()
for( id in allid ){
   #mapinfo of mismatch=0 
   M0 <- subset(map, N_mismatch==0 & ID==id)
   if( nrow(M0) == 0 ){ 
      print(paste0("WARNING! ", id," has no perfect match in the genome"))
      M0_map_info <- NA
   } else {
      M0_map_info <- paste(M0$map_info,collapse=":")
   }
   N_M0 <- nrow( subset(map, N_mismatch==0 & ID==id) )
   N_M1 <- nrow( subset(map, N_mismatch==1 & ID==id) )
   N_M2 <- nrow( subset(map, N_mismatch==2 & ID==id) )
   #dump <- data.frame(Primer_ID=id, M0_map_info, N_M0, N_M1, N_M2)
   dump <- data.frame(Primer_ID=id, N_M0, N_M1, N_M2)
      #2023-12-01: excluded M0_map_info since too long info sometimes
   res <- rbind(res,dump)
}

 #map info using 3-prime 15bp primer seq
map <- read.table(paste0(ODIR,"/g38_map.15.summary.txt"), header=T, sep="\t")

if (nrow(map) > 0) {
  map$map_info <- paste0(map$Chr, ":", map$Position, ":", map$Strand)
} else {
  # Create an empty data frame with the required structure
  map$map_info <- character(0)
}


map$map_info <- paste0(map$Chr, ":", map$Position, ":", map$Strand)
allid <- as.character(unique(map$ID))
res2 <- data.frame()
for( id in allid ){
   #mapinfo of mismatch=0 
   M0 <- subset(map, N_mismatch==0 & ID==id)
   if( nrow(M0) == 0 ){ 
      print(paste0("WARNING! ", id," has no perfect match in the genome"))
      M0_map_info_15bp <- NA
   } else {
      M0_map_info_15bp <- paste(M0$map_info,collapse=":")
   }
   N_M0_15bp <- nrow( subset(map, N_mismatch==0 & ID==id) )
   N_M1_15bp <- nrow( subset(map, N_mismatch==1 & ID==id) )
   N_M2_15bp <- nrow( subset(map, N_mismatch==2 & ID==id) )
   
   #dump <- data.frame(Primer_ID=id, M0_map_info_15bp, N_M0_15bp, N_M1_15bp, N_M2_15bp)
   dump <- data.frame(Primer_ID=id, N_M0_15bp, N_M1_15bp, N_M2_15bp)
      #2023-12-01: excluded M0_map_info_15bp since too long info
   res2 <- rbind(res2,dump)
}

 #merge and output
m <- merge(primer, res, by="Primer_ID")
m <- merge(m, res2, by="Primer_ID")

TB <- table(as.character(m$Name))
allnames <- names(TB[TB==2]) #the loci with valid in- and out-primers
out <- data.frame()
for(name in allnames){
   p_in <- subset(m, Name==name & DATA=="Inner")
   colnames(p_in)[-3] <- paste0(colnames(p_in)[-3],"_IN")
   p_out <- subset(m, Name==name & DATA=="Outer")
   colnames(p_out)[-3] <- paste0(colnames(p_out)[-3],"_OUT")
   dump <- merge(p_in,p_out,by="Name")
   out <- rbind(out, dump)
}
m <- merge(input ,out, by="Name")
OFILE=paste0(ODIR,"/target.primers.mapinfo_oneline.txt")
write.table(m, OFILE, col.names=T, row.names=F, append=F, quote=F, sep="\t")
