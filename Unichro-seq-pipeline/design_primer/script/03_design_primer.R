args <- commandArgs(trailingOnly = T)
allparameters <- as.character(args[1])
 # allparameters <- "Tm=62;Window=50;CVG=0.8;ODIR=analysis/2023-12-25/CVG0.8;IFILE=analysis/2023-12-25/CVG0.8/R11_v1.input"

library(openPrimeR)

allparameters <- unlist(strsplit(allparameters,";"))
param_names <- sapply(allparameters, function(x){
   x <- unlist(strsplit(x,"=")); return(x[1])
})
param <- sapply(allparameters, function(x){
   x <- unlist(strsplit(x,"=")); return(x[2])
})
names(param) <- param_names
print(param)

Tm <- as.numeric(param["Tm"])
   #Dist_snp <- as.numeric(param["Dist_snp"]) not present in this script
Window <- as.numeric(param["Window"])
CVG <-  as.numeric(param["CVG"])
ODIR <- as.character(param["ODIR"])
IFILE <- as.character(param["IFILE"])

print("Tm: melting temperature of primer")
   #print("Dist_snp: distance from snp to the primer starting point. minimum of allowed value = 1")
print("Window: window of site where we design primers")
print("CVG: 0-1. required.cvg in openPrimeR. Defines the desired coverage ratio of the templates. 1 indicates it design primers for all sites by relaxing the constrains. 0 indicates it design primers without relaxing the constrains.")
print("ODIR: the output directory")

 ###0, primer conditions
settings.xml <- system.file("extdata", "settings", "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
design.settings <- read_settings(settings.xml)
    #constraints: initial search condition
    #constraintLimits: conditions in the rescue search

constraints(design.settings)[["primer_length"]] <- c("min" = 18, "max" = 22) #same as default
constraintLimits(design.settings)[["primer_length"]] <- c("min" = 15, "max" = 25) #default 18-22

constraints(design.settings)[["melting_temp_range"]] <- c("min" = Tm - 3, "max" = Tm + 3) #default 55-65
constraintLimits(design.settings)[["melting_temp_range"]] <-  c("min" = Tm - 5, "max" = Tm + 5) #default 45-75

constraints(design.settings)[["melting_temp_diff"]] <- c("min" = 0, "max"=5) #default max = 5
constraintLimits(design.settings)[["melting_temp_diff"]] <- c("min" = 0, "max"=7.5) #default max = 7.5

constraints(design.settings)[["self_dimerization"]] <- c("min" = -15) #default -5 
constraintLimits(design.settings)[["self_dimerization"]] <- c("min" = -17) #default -7

constraints(design.settings)[["secondary_structure"]] <- c("min" = -11) #default -1 
constraintLimits(design.settings)[["secondary_structure"]] <- c("min" = -12) #default -2

constraints(design.settings)[["cross_dimerization"]] <- c("min" = -17) #default -7
constraintLimits(design.settings)[["cross_dimerization"]] <- c("min" = -19) #default -9

constraints(design.settings)[["primer_coverage"]] <- c("min"=1, "max"=1) #default min 1
constraintLimits(design.settings)[["primer_coverage"]] <- c("min"=1, "max"=1) #default min 1

conOptions(design.settings)[["allowed_mismatches"]] <- 0 #default 7
conOptions(design.settings)[["allowed_other_binding_ratio"]] <- 0 #default 1

cvg_constraints(design.settings) <- list() #stop coverage calculation


 ###1, Inner primer desing
 #read in input file
input_info <- read.table(IFILE,header=F,sep="\t")
if( ncol(input_info) != 4 ){  print("ERROR in the IFILE format!")  }
colnames(input_info) <- c("Varid","Peak_info","Strand_info","Dist_snp")

 #read in fasta file
fasta.file=paste0(ODIR,"/target.200bp.fa")
seq.df <- read_templates(fasta.file)

template.df.uni <- assign_binding_regions(
   seq.df,
   fw = c(1,50))
    #rev = c(1,50)
    #template.df.uni$Allowed_fw
    #the order is different from input_info

for( i in 1:nrow(template.df.uni)){
   Target_ID <- as.character(template.df.uni$Identifier[ i ])
   SNPID <- as.character(subset(template.df.uni, Identifier==Target_ID)$ID)
   SNPID <- gsub(">","",SNPID)
   Dist_snp <- subset(input_info, Varid==SNPID)$Dist_snp
   
   template.df.uni$Allowed_Start_fw[ i ] <- 200 - Dist_snp - Window + 1
   template.df.uni$Allowed_End_fw[ i ] <-  200 - Dist_snp
}

 # show(template.df.uni$Allowed_Start_fw)
 # show(template.df.uni$Allowed_End_fw)

idlist <- data.frame(
   Basic_Covered_Seqs=as.character(template.df.uni$Identifier),
   Name=gsub("^>","",template.df.uni$ID) )

 #search primer
optimal.primers <- design_primers(
   template.df.uni,
   required.cvg = CVG,
   mode.directionality = "fw",
   settings = design.settings)

inner_info <- optimal.primers$opti
OFILE=paste0(ODIR,"/inner_info.rds") #save all data before filtering
saveRDS(inner_info, file = OFILE)
inner_info <- inner_info[,c("Forward", "Basic_Covered_Seqs", "Identifier", "ID",
    "Basic_Binding_Position_Start_fw", "Basic_Binding_Position_End_fw",
   "primer_length_fw", "gc_clamp_fw", "gc_ratio_fw", "no_runs_fw", "no_repeats_fw",
   "Tm_C_fw","melting_temp_diff", "Self_Dimer_DeltaG","Cross_Dimer_DeltaG","Structure_deltaG",
   "primer_coverage", "Exclusion_Reason")]
inner_info <- merge(idlist, inner_info,by="Basic_Covered_Seqs")
inner_info$DATA <- "Inner"


 ###2, Outer primer desing
 #set primer sites
 #allowing 10 bp overlap between inner and outer primers
pass_inner_id <- inner_info$Basic_Covered_Seqs
template.df.uni2 <- template.df.uni[ template.df.uni$Identifier %in% pass_inner_id, ]
for( i in 1:nrow(template.df.uni2)){
   Target_ID <- as.character(template.df.uni2$Identifier[ i ])
   Inner_start_site <- as.numeric(subset(inner_info, Basic_Covered_Seqs==Target_ID)$Basic_Binding_Position_Start_fw)
   template.df.uni2$Allowed_Start_fw[ i ] <- Inner_start_site - Window + 10
   template.df.uni2$Allowed_End_fw[ i ] <- Inner_start_site + 10
}

 #widen acceptable Tm, only for rescue
 #constraintLimits(design.settings)[["melting_temp_range"]] <- c("min" = Tm - 10, "max" = Tm + 3)

 #search primer
optimal.primers <- design_primers(
   template.df.uni2,
   required.cvg = CVG,
   mode.directionality = "fw",
   settings = design.settings)

outer_info <- optimal.primers$opti

OFILE_outer=paste0(ODIR,"/outer_info.rds")
saveRDS(outer_info, file = OFILE_outer)

#This makes error unless you have already installed MELTING and  OligoArrayAux
outer_info <- outer_info[,c("Forward","Basic_Covered_Seqs", "Identifier", "ID",
    "Basic_Binding_Position_Start_fw", "Basic_Binding_Position_End_fw",
   "primer_length_fw", "gc_clamp_fw", "gc_ratio_fw", "no_runs_fw", "no_repeats_fw",
   "Tm_C_fw","melting_temp_diff", "Self_Dimer_DeltaG","Cross_Dimer_DeltaG","Structure_deltaG",
   "primer_coverage", "Exclusion_Reason")]
outer_info <- merge(idlist, outer_info,by="Basic_Covered_Seqs")
outer_info$DATA <- "Outer"

 #output
dump <- rbind(inner_info, outer_info)
c_order <- colnames(dump)
dump$Primer_ID <- paste0("P",1:nrow(dump))
dump <- dump[,c("Primer_ID", c_order)] #modify column order
OFILE=paste0(ODIR,"/target.primers")
write.table(dump, OFILE, col.names=T, row.names=F, append=F, quote=F, sep="\t")
