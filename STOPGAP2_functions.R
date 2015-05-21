#--------STOPGAP2 R functions ----------------#
# Author: Judong Shen

# Trim white spaces
trimWhiteSpace <- function (x) {
  sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

# Download GraspFullDataset2.zip data
# url <- c("https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip")
# file <- c("GraspFullDataset2.zip")
# download.file(url, file)

# grasp.import: use unix commands/R scripts to filter the grasp data
# Assume the input data is a *.zip file originally downloaded from the GRASP website.
grasp.import <- function(file = "GraspFullDataset2.zip", p.thres=1E-04)
{
  cat(paste("Start processing the ", file, " data ...", sep=""), "\n")
  system(paste("zcat", file,  "| head -n 1 > header.txt", sep=" "))
  vars <- names(read.delim("header.txt",
                           header = TRUE, comment.char = "", sep = "\t", as.is = TRUE))
  vars.keep <- c("NHLBIkey", "SNPid.dbSNP134.", "PMID", "LocationWithinPaper", "Pvalue",
                 "Phenotype", "DatePub", "Initial.Sample.Description", "Replication.Sample.Description",
                 "GWASancestryDescription", "TotalSamples.discovery.replication.", "TotalDiscoverySamples",
                 "Total.replication.samples")
  
  # Check whether all the variables we want to keep are in this version of GRASP data 
  stopifnot(all(vars.keep %in% vars))
  if (!all(vars.keep %in% vars)){
    cat("Not all variables we want to keep are in this version of GRASP data, please double check before running!", "\n")
  }
  
  no <- match(vars.keep, vars)
  no <- paste(no, collapse =",")
  p.no <- match("Pvalue", vars.keep)
  system(paste("zcat ", file, " | cut -f", no, " > tmp.txt", sep=""))
  system(paste("awk '$", p.no, "<=", p.thres,"' tmp.txt > grasp.txt", sep=""))
  system(paste("cut -f", no, " header.txt  > header1.txt", sep=""))
  system("cat header1.txt grasp.txt > grasp_analysis.txt")
  
  # Clean up
  unlink(c("header.txt","header1.txt", "tmp.txt", "grasp.txt"))
  
  # Read in the data
  grasp <- read.delim("grasp_analysis.txt",
                      na.strings = c("", "NA"), strip.white = TRUE)
  grasp <- grasp[-grep("^Gene expression", grasp$Phenotype),]
  grasp$PMID <- trimWhiteSpace(grasp$PMID)
  save(grasp, file = "grasp.RData", compress = TRUE)
  
  invisible(grasp)
}


# NHGRI.import: filter the NHGRI data
NHGRI.import <- function(file = "gwascatalog.txt", p.thres=1E-04)
{
  cat(paste("Start processing the NHGRI ", file, " data ...", sep=""), "\n")
  NHGRI <- read.delim(file,
                      header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)  
  NHGRI <- NHGRI[-grep("^[a-z | A-Z]", NHGRI$"p.Value"), ]
  NHGRI$"p.Value" <- as.numeric(NHGRI$"p.Value")
  NHGRI <- subset(NHGRI, p.Value <= p.thres)
  vars.keep <- c("PUBMEDID","SNPs", "p.Value", "Disease.Trait", "Initial.Sample.Description", "Replication.Sample.Description")
  NHGRI <- NHGRI[,vars.keep]
  names(NHGRI)[names(NHGRI)=="Initial.Sample.Description"] <- "Initial.Sample.Size"
  names(NHGRI)[names(NHGRI)=="Replication.Sample.Description"] <- "Replication.Sample.Size"
  
  # Remove the rows with haplotypes (multiple SNPs) instead of single SNPs
  NHGRI$SNPs <- trimWhiteSpace(NHGRI$SNPs)
  NHGRI$PUBMEDID <- trimWhiteSpace(NHGRI$PUBMEDID)
  nhgri <- NHGRI[-grep("[,|, ]", NHGRI$SNPs),]
  save(nhgri, file = "nhgri.RData", compress = TRUE)
  
  invisible(nhgri)
}

# gwasdb.import: filter the gwasdb data
gwasdb.import <- function(file = "gwasdb_20140812_snp_trait.gz", p.thres=1E-04)
{
  cat(paste("Start processing the gwasdb data ", file, "  ...", sep=""), "\n")
  gwasdb <- read.delim(gzfile(file),
                       header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  
  gwasdb <- gwasdb[-grep("^[a-z | A-Z]", gwasdb$P_VALUE), ]
  gwasdb$P_VALUE <- as.numeric(gwasdb$P_VALUE)  
  
  gwasdb <- subset(gwasdb, P_VALUE <= p.thres)
  gwasdb <- gwasdb[-grep("^Gene expression", gwasdb$GWAS_TRAIT),]
  vars.keep <- c("PMID","SNP_ID","ORI_SNP", "P_VALUE", "GWAS_TRAIT", "GWAS_INITIAL_SAMPLE_SIZE","SOURCE", "HPO_ID", "HPO_TERM", "DO_ID", "DO_TERM")
  gwasdb <- gwasdb[,vars.keep]
  gwasdb$SNP_ID <- trimWhiteSpace(gwasdb$SNP_ID)
  gwasdb$ORI_SNP <- trimWhiteSpace(gwasdb$ORI_SNP)
  gwasdb$PMID <- trimWhiteSpace(gwasdb$PMID)
  save(gwasdb, file = "gwasdb.RData", compress = TRUE)
  
  invisible(gwasdb)
}

# venn.plot: venn diagram plot
venn.plot <- function(id1, id2, id3){
  png(file="grasp_nhgri_gwasdb_venn.png",width=6.5,height=6.5,units="in",res=300)
  venn(list(GRASP=id1,NHGRI=id2,GWASDB=id3))
  dev.off()
}

# data.merge: merge the grasp, nhgri and gwasdb datasets by first using nhgri records, supplementing 
#              by grasp and gwasdb datasets respectively
# Merge order: NHGRI > GRASP > GWASDB
data.merge <- function(nhgri, grasp, gwasdb){
  cat("Start merging the data from three sources: nhgri, grasp, and gwasdb ...", "\n")
  
  nhgri <- rename(nhgri, c(SNPs = 'snp_id', 
                           p.Value = "pvalue",
                           Disease.Trait = 'disease',
                           Initial.Sample.Size = 'Initial.Sample',
                           Replication.Sample.Size = "Replication.Sample"))
  nhgri$pvalue <- as.numeric(format(nhgri$pvalue, digits=3, scientific=T))  
  nhgri$Source <- "nhgri"
  
  # Remove redundant rows in nhgri data ("PUBMEDID", "snp_id", "pvalue", "disease", "Initial.Sample", "Replication.Sample" "Source" )
  nhgri <- nhgri[!duplicated(nhgri),]
  
  grasp <- rename(grasp, c(SNPid.dbSNP134. = 'snp_id', 
                           PMID = 'PUBMEDID',
                           Pvalue = "pvalue",
                           Phenotype = 'disease',
                           Initial.Sample.Description = 'Initial.Sample',
                           Replication.Sample.Description = 'Replication.Sample'))
  grasp$snp_id <- paste("rs",grasp$snp_id,sep="")
  grasp$Source <- "grasp"
  grasp <- subset(grasp, !is.na(pvalue))
  grasp$pvalue <- as.numeric(format(grasp$pvalue, digits=3, scientific=T)) 
  
  gwasdb <- rename(gwasdb, c(PMID = 'PUBMEDID',
                             SNP_ID = "gwasdb_SNP_ID",
                             ORI_SNP = 'snp_id', 
                             P_VALUE = "pvalue",
                             GWAS_TRAIT = 'disease',
                             GWAS_INITIAL_SAMPLE_SIZE = 'Initial.Sample',
                             SOURCE = "Source_gwasdb"))

  # Maually change some of the mis-spelled pvalue here
  gwasdb$pvalue[gwasdb$pvalue=="1.49 E-9"] <- "1.49E-9"
  gwasdb$pvalue[gwasdb$pvalue=="1.05 E-14"] <- "1.05E-14"
  gwasdb$pvalue[gwasdb$pvalue=="1.3E?\\05"] <- "1.3E-05"
  gwasdb$pvalue[gwasdb$pvalue=="1.78E?05"] <- "1.78E-05"
  gwasdb <- subset(gwasdb, pvalue!="")
  gwasdb$pvalue <- as.numeric(format(as.numeric(gwasdb$pvalue), digits=3, scientific=T)) 
  
  grasp$id <- paste(grasp$snp_id, grasp$PUBMEDID,sep="_")
  nhgri$id <- paste(nhgri$snp_id, nhgri$PUBMEDID,sep="_")
  gwasdb$id <- paste(gwasdb$snp_id, gwasdb$PUBMEDID,sep="_")
  id1 <- unique(grasp$id); id2=unique(nhgri$id); id3=unique(gwasdb$id)
  
  venn.plot(id1=id1, id2=id2, id3=id3)
  
  # Subset of grasp data with id (SNP_PMID) not in nhgri data
  grasp1 <- subset(grasp, !(id %in% nhgri$id))
  grasp1$Source <- "grasp"
  # Remove redundant rows in grasp1 data 
  grasp1 <- grasp1[!duplicated(grasp1),]
  
  # Subset of grasp data with id (SNP_PMID) not in nhgri and grasp
  gwasdb1 <- subset(gwasdb, !(id %in% c(nhgri$id, grasp$id)))
  gwasdb1$Source <- "gwasdb"
  # Remove redundant rows in grasp1 data 
  gwasdb1 <- gwasdb1[!duplicated(gwasdb1),]  
  
  data <- merge(nhgri, grasp1, 
                by = c("snp_id","PUBMEDID", "pvalue","disease","Initial.Sample",
                       "Replication.Sample","Source"), all.x = TRUE, all.y = TRUE)
  data <- merge(data, gwasdb1,
                by = c("snp_id","PUBMEDID", "pvalue","disease","Initial.Sample","Source"),
                all.x = TRUE, all.y = TRUE)
  data$snp_id <- trimWhiteSpace(data$snp_id)
  data <- data[,-match(c("id.x","id.y", "id"), names(data))]
  
  #write.table(data,
  #            "stopgap_grasp_nhgri_gwasdb_merged.txt",
  #            sep = "\t", row.names = FALSE, quote = F, na = "")
  save(data, file = "stopgap_3sources.RData", compress = TRUE)
  
  invisible(data)
  
}


# gwas.gap: identify the additional GWAS results compared with the previous version:
gwas.gap <- function(new.gwas="stopgap_3sources.RData", old.gwas="../../STOPGAP_pipeline/Data/stopgap_3sources.RData"){
  load(new.gwas)
  new.gwas <- data
  load(old.gwas)
  old.gwas <- data
  if (all(old.gwas$snp_id %in% new.gwas$snp_id)){
    gwas.gap <- subset(new.gwas, !(snp_id %in% old.gwas$snp_id))
    write.table(gwas.gap,
               "gwas.gap.txt",
               sep = "\t", row.names = FALSE, quote = F, na = "")
    save(gwas.gap, file = "gwas.gap.RData", compress = TRUE)    
  } else {
    cat("Not all snp_id in the previous version of GWAS data are in this version, please double check before running!", "\n")
  }

  invisible(gwas.gap)
}

# gwas.filter: filter the gwas data based on disease names and PUBMEDID
gwas.filter <- function(data){
  # Remove the gwas data with some disease names that will not be used 
  no <- grep("methylation levels", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("differential exon", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("differential splicing", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("transcript initiation", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("transcript termination", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("gene expression", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("methylation qtl", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("qtl", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("exon skipping", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("differential", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("desialylated glycan peak", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("complex transcript isoform variation", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("recombination", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("intron retention", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("transcript levels", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("allele-specific methylation", data$disease)
  if (length(no)>0){
    data <- data[-no,]
  }   
  
  # Further remove records in some publications
  data <- subset(data, !(PUBMEDID %in% c("21886157",  # Serum ratio of ...
                                         "18403759",  # YKL-40 levels
                                         "21572414")))  # "Urinary metabolites"
  #   write.table(data,
  #               "stopgap_grasp_nhgri_gwasdb_merged_filtered.txt",
  #               sep = "\t", row.names = FALSE, quote = F, na = "")
  #   save(data, file = "stopgap_3sources_filtered.RData", compress = TRUE)
  data
  
}

# new.gwassnps: identify the unique new GWAS SNPs compared with the previous version:
gwas.snps <- function(gwas.file="stopgap_3sources.RData"){
  load(gwas.file)
  snps <- sort(unique(data$snp_id))
  writeLines(snps, "stopgap2_SNPs.txt", sep = "\n", useBytes = FALSE)
  
  invisible(snps)
}

# Coordinate lookup and LD calculation based on the 1KG data for the new SNPs identified in this version
# run.ld: Calculate the LD SNPs based on the 
# rslist is the list of rs IDs with # prefixed comment lines allowed
run.ld <- function(){
  # Input SNPs: stopgap2_SNPs.txt
  
  system("mkdir STOPGAP2_LDResults")
  
  #(1) Coordinate lookup
  system("./STOPGAP2_LD/Get_rsID_coord.py -o ./STOPGAP2_LDResults -r ./stopgap2_SNPs.txt") 
    
  #(2) LD Calculation
  system("./STOPGAP2_LD/LDabq.sh ./STOPGAP2_LDResults ./STOPGAP2_LDResults")
  
}

# Update rsIDs in the gwas data (data) to dbSNP141 version by using the ./STOPGAP2_LDResults/rsID_Coordinates.txt file
# Write out the new GWAS data as "stopgap_3sources_dbSNP141.RData"
rsID.update <- function(gwas.file="stopgap_3sources.RData", 
                        coor.file="./STOPGAP2_LDResults/rsID_Coordinates.txt")
                        #"T:/statgen4/STOPGAP_LD/v2.1/rsID_Coordinates_supplemented.txt"
{
  load(gwas.file)
  gwas.data <- gwas.filter(data)
  gwas.data$disease <- tolower(trimWhiteSpace(gwas.data$disease))
  
  rsid.coor <- read.delim(coor.file,
                         header = TRUE, comment.char = "", sep = "\t", as.is = TRUE) 
  names(rsid.coor)[names(rsid.coor)=="X.rsID"] <- "rsID "
  gwas.data$rsid <- rsid.coor$rsID[match(gwas.data$snp_id, rsid.coor$inputrsIDs )]
  names <- names(gwas.data)[!(names(gwas.data) %in% c("snp_id", "rsid"))]
  gwas.data <- gwas.data[,c("snp_id", "rsid", names)]
  
  save(gwas.data, file = "stopgap_3sources_dbSNP141.RData", compress = TRUE) # 176159 unique GWAS SNPs
  
  invisible(gwas.data)
}


# ChIA.PET.import: Process ChIA-PET files and merge them together
ChIA.PET.import <- function(ChIA.PET.path="./ChIA-PET", 
                            gencode.file="gencode_v20.RData",
                            win.size=1000){
  cat("Start processing the ChIA-PET data ...", "\n")
  files <- list.files(ChIA.PET.path, pattern="*.bed$")
  files <- files[grep("^Pol2", files)]
  n <- length(files)
  
  load(gencode.file)
  gencode$Gene <- as.character(gencode$Gene)
  out <- list()
  for (i in 1:n){
    file <- files[i]
    cat(file, "\n")
    file1 <- paste(ChIA.PET.path, file, sep="/")
    ChIA.data <- read.delim(file1,
                            header = FALSE, comment.char = "", sep = "\t", as.is = TRUE)
    names(ChIA.data) <- c("chr1", "startpos1", "endpos1", "region")
    ChIA.data$chr1 <- gsub("chr", "", ChIA.data$chr)
    r.chr <- unlist(lapply( strsplit(ChIA.data$region,":"),function(x)x[[1]]))
    ChIA.data$chr2 <- gsub("chr", "", r.chr)
    r.range <- unlist(lapply( strsplit(ChIA.data$region,":"),function(x)x[[2]]))
    ChIA.data$startpos2 <- unlist(lapply( strsplit(r.range,"-"),function(x)x[[1]]))
    ChIA.data$endpos2 <- unlist(lapply( strsplit(r.range,"-"),function(x)x[[2]]))
    ChIA.data$region <- NULL
    
    P.Genes <- unlist(apply(ChIA.data, 1, function(x){
      gencode.sub1 <- subset(gencode, chrom==x[1] & start > (as.numeric(x[2]) - win.size) & start < (as.numeric(x[3]) + win.size) )
      p.genes <- gencode.sub1$Gene
      n.p.genes <- length(p.genes)
      pgenes <- NA
      if (n.p.genes>0){
        pgenes <- paste(p.genes, collapse="__")
      } 
      pgenes
    }))
    
    D.Genes <- unlist(apply(ChIA.data, 1, function(x){
      # x <- as.character(ChIA.data[7,])
      gencode.sub2 <- subset(gencode, chrom==x[4] & start > (as.numeric(x[5]) - win.size) & start < (as.numeric(x[6]) + win.size) )
      d.genes <- gencode.sub2$Gene
      n.d.genes <- length(d.genes)
      dgenes <- NA
      if (n.d.genes>0){
        dgenes <- paste(d.genes, collapse="__")
      }   
      dgenes
    }))
    
    ChIA.data1 <- ChIA.data2 <- ChIA.data
    ChIA.data1$pgenes <- P.Genes
    ChIA.data2$dgenes <- D.Genes
    ChIA.data1 <- subset(ChIA.data1, !is.na(pgenes))
    ChIA.data2 <- subset(ChIA.data2, !is.na(dgenes))
    ChIA.data1 <- ChIA.data1[,c("chr2", "startpos2", "endpos2", "chr1", "startpos1", "endpos1", "pgenes")]
    names <- c("dChr", "dStart", "dEnd", "pChr", "pStart", "pEnd", "gene")
    names(ChIA.data1) <- names
    #   Divide out rows with more than one gene
    genes <- strsplit(as.character(ChIA.data1$gene), "__")
    ChIA.data1 <- data.frame(ChIA.data1[rep(1:nrow(ChIA.data1), unlist(lapply(genes, length))), -ncol(ChIA.data1)], 
                             genes = factor(unlist(genes)))
    names(ChIA.data2) <- names
    genes <- strsplit(as.character(ChIA.data2$gene), "__")
    ChIA.data2 <- data.frame(ChIA.data2[rep(1:nrow(ChIA.data2), unlist(lapply(genes, length))), -ncol(ChIA.data2)], 
                             genes = factor(unlist(genes)))
    ChIA.data <- rbind(ChIA.data1, ChIA.data2)
    ChIA.data$dChr <- factor(ChIA.data$dChr, levels=c(as.character(1:22), "X"))
    ChIA.data$pChr <- factor(ChIA.data$pChr, levels=c(as.character(1:22), "X"))
    ChIA.data$dStart <- as.numeric(ChIA.data$dStart)
    ChIA.data$dEnd <- as.numeric(ChIA.data$dEnd)
    ChIA.data$pStart <- as.numeric(ChIA.data$pStart)
    ChIA.data$pEnd <- as.numeric(ChIA.data$pEnd)
    ChIA.data <- ChIA.data[order(ChIA.data$dChr, ChIA.data$dStart, ChIA.data$dEnd, ChIA.data$pChr, ChIA.data$pStart, ChIA.data$pEnd, ChIA.data$gene),]
    
    ChIA.data$CellType <- file
    #write.table(ChIA.data, paste(file,".txt",sep=""),
    #            sep = "\t", row.names = FALSE, quote = F, na = "")    
    
    out[[i]] <- ChIA.data
  }
  out <- do.call("rbind", out)
  
  CHIA.PET <- out
  #write.table(CHIA.PET, "CHIA_PET_Pol2.txt",
  #            sep = "\t", row.names = FALSE, quote = F, na = "") 
  system("mkdir CHIA_PET")
  save(CHIA.PET, file = "./CHIA_PET/CHIA.PET.RData", compress = TRUE)  
  #rm(out)
  
  invisible(CHIA.PET)
}

# ld.plusTarget: include the target GWAS rsIDs as part of the LD SNPs as well
ld.plusTarget <- function(ld.data){
  if ("R.2" %in% names(ld.data)){
    names(ld.data)[names(ld.data)=="R.2"] <- "r2"
  }
  ld.data <- ld.data[,c("CHR1", "POS1", "ID1", "REF1", "ALT1", "CHR2", "POS2", "ID2", "REF2", "ALT2", "r2" )]
  ld.data1 <- ld.data[,c("CHR1", "POS1", "ID1", "REF1", "ALT1")]
  ld.data1 <- ld.data1[!duplicated(ld.data1),]
  ld.data2 <- ld.data1
  names(ld.data2) <- c("CHR2", "POS2", "ID2", "REF2", "ALT2")
  ld.data2$r2 <- 1
  ld.data.target <- cbind(ld.data1, ld.data2)
  ld.data <- rbind(ld.data, ld.data.target)
  ld.data <- ld.data[order(ld.data$CHR1, ld.data$POS1, ld.data$ID1, ld.data$CHR2, ld.data$POS2, ld.data$ID2),]
  names(ld.data) <- c("CHR.gwas", "POS.gwas", "SNP.gwas", "REF.gwas", "ALT.gwas", "CHR.ld", "POS.ld", "SNP.ld", "REF.ld", "ALT.ld", "r2" )
  ld.data
}

# Find all variants in LD with GWAS variants at r2 > 0.5 - GWAS, previous version of LD data, the path to the new additional LD data
ld.snps <- function(gwas.data="stopgap_3sources_dbSNP141.RData",
                    ld.path="./STOPGAP2_LDResults"){
  cat("Start finding all the SNPs with LD r2 > 0.5 ...", "\n")
  load(gwas.data)
  
  gwas.ld.snps <- gwas.data[,c("snp_id", "rsid")]
  gwas.ld.snps <-  gwas.ld.snps[!duplicated(gwas.ld.snps),]
  
  ld.files <- list.files(ld.path, pattern="*.labeled.ld$") 
  n <- length(ld.files)
  chrs <- unlist(lapply( strsplit(ld.files,"\\."),function(x)x[[1]]))
  chrs <- gsub("chr", "", chrs)
  m.no <- match(c(as.character(1:22), "X"), chrs)  
  
  ld.snps.r2 <- list() 
  k <- 0
  for (i in m.no){
    k <- k + 1
    cat(k, "\n")
    file <- ld.files[i]
    file <- paste(ld.path, file, sep="/")
    ld.data <- read.delim(file,
                          header = TRUE, comment.char = "", sep = "", as.is = TRUE)
    names(ld.data)[names(ld.data)=="R.2"] <- "r2"
    ld.data <- ld.plusTarget(ld.data)
    
    if (nrow(ld.data)>0){
      # ld data includes more SNPs than in the gwas data, so filter down it first. 
      ld.data <- subset(ld.data, SNP.gwas %in% gwas.ld.snps$rsid)
    } else {
      ld.data <- NA
    }
    
    ld.snps.r2[[k]] <- ld.data
  }
  ld.snps.r2 <- do.call("rbind", ld.snps.r2) 
  system("mkdir ./STOPGAP2_LDResults/gwas_LD")
  save(ld.snps.r2, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.RData", compress = TRUE)
  
  ld.snps <- ld.snps.r2[,c("CHR.ld", "POS.ld", "SNP.ld")]
  ld.snps <- ld.snps[!duplicated(ld.snps),]
  save(ld.snps, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.RData", compress = TRUE)  
  
  invisible(ld.snps.r2)
}

# updata.ld.data: updata ld data:
# Remove SNP.ld with two positions in both ld.snps and ld.snps.r2
# Use chr:pos to replace the missing ld SNP ids in both ld.snps and ld.snps.r2
updata.ld.data <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.RData", 
                           ld.r2.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.RData"){
  # "T:/statgen4/STOPGAP_LD/v4/LDcalc/gwas_LD/GWAS_LD_SNPs.RData
  load(ld.file)
  load(ld.r2.file)
  
  # Remove SNP.ld with two positions in both ld.snps and ld.snps.r2
  snps.rm <- names(table(ld.snps$SNP.ld))[which(table(ld.snps$SNP.ld)==2)]
  ld.snps.r2 <- subset(ld.snps.r2, !(SNP.ld %in% snps.rm))
  ld.snps <- subset(ld.snps, !(SNP.ld %in% snps.rm))
  
  # Use chr:pos to replace the missing ld SNP ids in both ld_snps and ld.snps.r2
  miss.no <- which(ld.snps$SNP.ld==".")
  ld.snps$SNP.ld[miss.no] <- paste(ld.snps$CHR.ld[miss.no], ld.snps$POS.ld[miss.no], sep=":")  
  miss.no1 <- which(ld.snps.r2$SNP.ld==".")
  ld.snps.r2$SNP.ld[miss.no1] <- paste(ld.snps.r2$CHR.ld[miss.no1], ld.snps.r2$POS.ld[miss.no1], sep=":")  
  save(ld.snps.r2, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.RData", compress = TRUE)
  save(ld.snps, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.RData", compress = TRUE)  
  invisible(ld.snps.r2)
}

# Update the ld.snps.r2 and ld_snps data by including the GWAS SNPs data which don't have any LD information there
update.ld.snps.r2 <- function(gwas.file="stopgap_3sources_dbSNP141.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.RData",
                              rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt"){
  cat("Update the ld.snps.r2 data by including the (154526-147706) GWAS SNPs data there ...", "\n")
  load(gwas.file)
  load(ld.file)
  coord <- read.delim(rsid.file,
                      header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  names(coord)[names(coord)=="X.rsID"] <- "rsID"
  coord$chr <- unlist(lapply( strsplit(coord$GRCh37_Coordinate,":"),function(x)x[[1]]))
  coord$pos <- unlist(lapply( strsplit(coord$GRCh37_Coordinate,":"),function(x)x[[2]]))
  coord <- coord[!duplicated(coord),]
  
  # gwas data with rsid not in the ld.snps.r2 data
  gwas.data1 <- subset(gwas.data, !(rsid %in% unique(ld.snps.r2$SNP.ld)) )
  gwas.data1$chr <- coord$chr[match(gwas.data1$rsid, coord$rsID)]
  gwas.data1$pos <- coord$pos[match(gwas.data1$rsid, coord$rsID)]
  gwas.data1$REF <- coord$REF[match(gwas.data1$rsid, coord$rsID)]
  gwas.data1$ALT <- coord$ALT[match(gwas.data1$rsid, coord$rsID)]
  gwas.data2 <- gwas.data1[,c("chr", "pos", "rsid", "REF", "ALT")]
  gwas.data2 <- cbind(gwas.data2, gwas.data2)
  gwas.data2$r2 <- 1
  names(gwas.data2) <- names(ld.snps.r2)
  ld.snps.r2 <- rbind(ld.snps.r2, gwas.data2)
  save(ld.snps.r2, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.RData", compress = TRUE)
  
  ld.snps <- ld.snps.r2[,c("CHR.ld", "POS.ld", "SNP.ld")]
  ld.snps <- ld.snps[!duplicated(ld.snps),]
  save(ld.snps, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.RData", compress = TRUE)  
  
  invisible(ld.snps.r2)
}

# Update the ld.snps.r2 and ld_snps data by removing the LD SNPs more than 500kb away from the GWAS SNPs,
# If any GWAS SNPs get deleted, put them back by including r2=1 there
update1.ld.snps.r2 <- function(gwas.file="stopgap_3sources_dbSNP141.RData",
                              ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.RData",
                              dis = 500*1000){
  cat("Update the ld.snps.r2 data by removing the LD SNPs more than 500kb away from the GWAS SNPs ...", "\n")
  load(gwas.file)
  load(ld.file)
  ld.snps.r2$POS.gwas <- as.numeric(ld.snps.r2$POS.gwas)
  ld.snps.r2$POS.ld <- as.numeric(ld.snps.r2$POS.ld)
  ld.snps.r2 <-  subset(ld.snps.r2, abs(POS.gwas-POS.ld)<=dis)
  ld.snps.r2 <- ld.snps.r2[!duplicated(ld.snps.r2),]
  save(ld.snps.r2, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.RData", compress = TRUE)
  
  ld.snps <- ld.snps.r2[,c("CHR.ld", "POS.ld", "SNP.ld")]
  ld.snps <- ld.snps[!duplicated(ld.snps),]
  save(ld.snps, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", compress = TRUE)  
  
  gwas.r2 <- merge(gwas.data[,c("rsid", "PUBMEDID", "pvalue")], ld.snps.r2[,c("SNP.gwas", "SNP.ld", "r2", "CHR.ld", "POS.ld")], 
                   by.x = "rsid", by.y="SNP.gwas", all = TRUE)
  names(gwas.r2) <- c("snp.gwas", "pubmedid", "pvalue", "snp.ld", "r2", "chr.ld", "pos.ld")
  gwas.r2 <- gwas.r2[!duplicated(gwas.r2),]
  system("mkdir LocusClustering")
  save(gwas.r2, file = "./LocusClustering/gwas_r2.RData", compress = TRUE)  
  
  invisible(ld.snps.r2)
}

# Update the ld.snps.r2 data by adding the AF information from 1KG phase I data ...
r2.frq <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.RData",
                   frq.file="./1KG_AF/gws.frq.RData"){
  cat("Update the ld.snps.r2 data by adding the AF information from 1KG phase I data ...", "\n")
  load(ld.file)
  load(frq.file)
  
  # frq.1kg
  id.no <- match("ID", names(frq.1kg))
  ld.snps.r2 <- merge(ld.snps.r2, frq.1kg[,-id.no], 
                      by.x = c("CHR.ld","POS.ld", "REF.ld","ALT.ld"),
                      by.y=c("CHROM", "POS", "REF", "ALT"), all.x = TRUE, all.y = FALSE)
  
  save(ld.snps.r2, file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData", compress = TRUE)
  
  invisible(ld.snps.r2)
}

# ld.rdb.dhscor: add the Cat.rdb to the ld snps &
# For each LD SNP, (1) identify each dhscor where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.rdb.dhscor <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                      rdb.file="./ENCODE/regulomedb14.RData",
                      dhscor.file="./ENCODE/dhscor.RData",
                      mpos.rng.p = c("Chr.pdhs", "Start.pdhs", "Stop.pdhs"),
                      mpos.rng.d = c("Chr.ddhs", "Start.ddhs", "Stop.ddhs")){
  load(ld.file)
  load(rdb.file)
  
  # add the Cat.rdb to the ld snps 
  ld.snps$Cat.rdb <- rdb$Cat.rdb[match(ld.snps$SNP.ld, rdb$SNP)]  
  
  load(dhscor.file)
  dhscor$Chr.ddhs <- gsub("chr", "", dhscor$Chr.ddhs)
  dhscor$Chr.pdhs <- gsub("chr", "", dhscor$Chr.pdhs)
  dhscor$Gene <- as.character(dhscor$Gene)
  dhscor.p <- dhscor[,c("Chr.pdhs", "Start.pdhs", "Stop.pdhs", "Gene")]
  dhscor.d <- dhscor[,c("Chr.ddhs", "Start.ddhs", "Stop.ddhs", "Gene", "Cor")]
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    # promoter
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    dhscor.p.i <- subset(dhscor.p, Chr.pdhs==i)
    dhscor.p.i$Cor <- 1
    
    ## Expand range as one row per position for dhscor.p.i
    dhscor.p.i$Length <- dhscor.p.i[, mpos.rng.p[3]] - dhscor.p.i[, mpos.rng.p[2]] + 1
    dhscor.p.i.exp <- data.frame(Position = as.vector(unlist(apply(dhscor.p.i[, mpos.rng.p[2:3]], 1,
                                                                   function(x) x[1]:x[2]))),
                                 Gene = rep(as.character(dhscor.p.i$Gene), dhscor.p.i$Length),
                                 Cor = rep(dhscor.p.i$Cor, dhscor.p.i$Length))
    dhscor.p.i.exp <- dhscor.p.i.exp[!duplicated(dhscor.p.i.exp),]
    tmp <- merge(ld.snps.i, dhscor.p.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp1 <- tmp[!duplicated(tmp),]
    tmp1$type <- NA
    tmp1$type[!is.na(tmp1$Cor)] <- "p"
    
    # distal
    dhscor.d.i <- subset(dhscor.d, Chr.ddhs==i)
    
    ## Expand range as one row per position for dhscor.d.i
    dhscor.d.i$Length <- dhscor.d.i[, mpos.rng.d[3]] - dhscor.d.i[, mpos.rng.d[2]] + 1
    dhscor.d.i.exp <- data.frame(Position = as.vector(unlist(apply(dhscor.d.i[, mpos.rng.d[2:3]], 1,
                                                                   function(x) x[1]:x[2]))),
                                 Gene = rep(as.character(dhscor.d.i$Gene), dhscor.d.i$Length),
                                 Cor = rep(dhscor.d.i$Cor, dhscor.d.i$Length))
    dhscor.d.i.exp <- dhscor.d.i.exp[!duplicated(dhscor.d.i.exp),]
    tmp <- merge(ld.snps.i, dhscor.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$Cor)] <- "d"
    
    tmp <- merge(tmp1, tmp2, 
                 by = c("CHR.ld","POS.ld", "SNP.ld", "Cat.rdb","Gene",  "Cor", "type"), 
                 all.x = TRUE, all.y = TRUE)
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  gwas.ld.rdb.dhscor <- res
  system("mkdir VarGeneMapping")
  save(gwas.ld.rdb.dhscor, file = "./VarGeneMapping/gwas.ld.rdb.dhscor.RData", compress = TRUE) 
  
  invisible(res)  
  
}

# ----- Merge the CHIA-PET data
# ld.chiapet: For each LD SNP, (1) indentify each chiapet where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.chiapet <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                       chiapet.file="./CHIA_PET/CHIA.PET.RData",
                       mpos.rng.p = c("pChr", "pStart", "pEnd"),
                       mpos.rng.d = c("dChr", "dStart", "dEnd")){
  cat("Start merging CHIA-PET data...", "\n")
  load(ld.file)
  load(chiapet.file)
  names(CHIA.PET)[names(CHIA.PET)=="genes"] <- "Gene"
  CHIA.PET$CellType <- gsub(".bed", "", CHIA.PET$CellType)
  CHIA.PET$Gene <- as.character(CHIA.PET$Gene)
  CHIA.PET.p <- CHIA.PET[,c("pChr", "pStart", "pEnd", "Gene", "CellType")]
  CHIA.PET.d <- CHIA.PET[,c("dChr", "dStart", "dEnd", "Gene", "CellType")]
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    # promoter
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    CHIA.PET.p.i <- subset(CHIA.PET.p, pChr==i)
    
    ## Expand range as one row per position for chiapet.p.i
    CHIA.PET.p.i$Length <- CHIA.PET.p.i[, mpos.rng.p[3]] - CHIA.PET.p.i[, mpos.rng.p[2]] + 1
    CHIA.PET.p.i.exp <- data.frame(Position = as.vector(unlist(apply(CHIA.PET.p.i[, mpos.rng.p[2:3]], 1,
                                                                     function(x) x[1]:x[2]))),
                                   Gene = rep(as.character(CHIA.PET.p.i$Gene), CHIA.PET.p.i$Length),
                                   CellType = rep(CHIA.PET.p.i$CellType, CHIA.PET.p.i$Length))
    CHIA.PET.p.i.exp <- CHIA.PET.p.i.exp[!duplicated(CHIA.PET.p.i.exp),]
    tmp <- merge(ld.snps.i, CHIA.PET.p.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp1 <- tmp[!duplicated(tmp),]
    tmp1$type <- NA
    tmp1$type[!is.na(tmp1$CellType)] <- "p"
    
    # distal
    CHIA.PET.d.i <- subset(CHIA.PET.d, dChr==i)
    
    ## Expand range as one row per position for chiapet.d.i
    CHIA.PET.d.i$Length <- CHIA.PET.d.i[, mpos.rng.d[3]] - CHIA.PET.d.i[, mpos.rng.d[2]] + 1
    CHIA.PET.d.i.exp <- data.frame(Position = as.vector(unlist(apply(CHIA.PET.d.i[, mpos.rng.d[2:3]], 1,
                                                                     function(x) x[1]:x[2]))),
                                   Gene = rep(as.character(CHIA.PET.d.i$Gene), CHIA.PET.d.i$Length),
                                   CellType = rep(CHIA.PET.d.i$CellType, CHIA.PET.d.i$Length))
    CHIA.PET.d.i.exp <- CHIA.PET.d.i.exp[!duplicated(CHIA.PET.d.i.exp),]
    tmp <- merge(ld.snps.i, CHIA.PET.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$CellType)] <- "d"
    
    # Remove rows with all NAs for the "Gene", "CellType" and  "type" columns before merging
    #tmp2$id <- paste(tmp2$Gene, tmp2$CellType, tmp2$type, sep="_")
    #tmp3 <- subset(tmp2, id!="NA_NA_NA")
    
    tmp <- merge(tmp1, tmp2, 
                 by = c("CHR.ld","POS.ld", "SNP.ld", "Gene",  "CellType", "type"), 
                 all.x = TRUE, all.y = TRUE)
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  gwas.ld.chiapet <- res
  save(gwas.ld.chiapet, file = "./VarGeneMapping/gwas.ld.chiapet.RData", compress = TRUE) 
  
  invisible(res)  
  
}

#eQTL.GTEx: process the GTEX.data.RData
eQTL.GTEx <- function(eQTL.GTEx.file="./eQTL/GTEx/GTEx.data.RData",
                       eQTL.GTEx.sz="./eQTL/GTEx/GETx_Tissue_SampleSize.txt"){
  cat("Start processing the GTEX.data.RData...", "\n")
  load(eQTL.GTEx.file)  
  GTEx.sz <- read.delim(eQTL.GTEx.sz,
                      na.strings = c("", "NA"), strip.white = TRUE)
  GTEx.data$SampleSize <- GTEx.sz$SampleSize[match(GTEx.data$tissue, GTEx.sz$Tissue)]
  
  eQTL <- GTEx.data[,c( "snp", "symbol", "tissue",  "SampleSize", "p.value", "fdr",
                        "nom_thresh", "min.p.", "empp", "ks","n")]
  names(eQTL) <- c("SNP","GeneSymbol", "Tissue", "SampleSize", "P.value", "FDR", "NomThresh", "Min.p", "EmpP", "KS","N")
  names(eQTL) <- paste(names(eQTL), "GTEx", sep=".")
  names(eQTL)[names(eQTL)=="SNP.GTEx"] <- "SNP"
  #eQTL$eQTLSource <- "GTEx"
  
  return(eQTL)
  
}

# ld.eQTL: For each LD SNP, (1) indentify each eQTL where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.eQTL <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                    eQTL.file="./eQTL/eQTL_UChichago/eQTL_UChicago_STOPGAP.RData",
                    eQTL.GTEx.file="./eQTL/GTEx/GTEx.data.RData",
                    eQTL.GTEx.sz="./eQTL/GTEx/GETx_Tissue_SampleSize.txt"){
  cat("Start merging eQTL data...", "\n")
  load(ld.file)
  load(eQTL.file)
  eQTL <- eQTL.uc2[,c( "Ref1", "SNP", "RegionGene", "Score", "ScoreType", "Tissue", "SampleSize", "Sample")]
  # The RegionGene column actually includes all genes
  # Change the "-" in the SNP column to ":" since it is ":" in the ld.snps data
  eQTL$SNP <- as.character(eQTL$SNP)
  m.no <- grep("-", eQTL$SNP)
  eQTL$SNP[m.no] <- gsub("-", ":", eQTL$SNP[m.no])
  #eQTL$eQTLSource <- "UChicago"
  
  eQTL <- subset(eQTL, SNP %in% ld.snps$SNP.ld)
  names(eQTL) <- paste(names(eQTL),"UCh",sep=".")
  names(eQTL)[names(eQTL)=="SNP.UCh"] <- "SNP"
  
  eQTL.gtex <- eQTL.GTEx(eQTL.GTEx.file="./eQTL/GTEx/GTEx.data.RData",
                        eQTL.GTEx.sz="./eQTL/GTEx/GETx_Tissue_SampleSize.txt")  
  eQTL.gtex <- subset(eQTL.gtex, SNP %in% ld.snps$SNP.ld)
  
  # For eQTL and eQTL.gtex, only keep the "SNP", "Gene", "Ref", "Tissue", "ScoreType", "Score" columns
  # And then combine them together
  eQTL.uc <- eQTL[,c("SNP", "RegionGene.UCh", "Ref1.UCh", "Tissue.UCh", "ScoreType.UCh", "Score.UCh")]
  names(eQTL.uc) <- c("SNP", "Gene", "Ref", "Tissue", "ScoreType", "Score")
  eQTL.uc$Ref <- paste("Uchicago", eQTL.uc$Ref, sep=".")
  
  eQTL.gtex1 <- eQTL.gtex[,c("SNP", "GeneSymbol.GTEx", "Tissue.GTEx", "P.value.GTEx")]
  eQTL.gtex1$ScoreType <- "-log10(P)"
  eQTL.gtex1$Score <- -log10(eQTL.gtex1$P.value.GTEx)
  eQTL.gtex1$Ref <- "GTEx"
  eQTL.gtex1 <- eQTL.gtex1[,c("SNP", "GeneSymbol.GTEx", "Ref", "Tissue.GTEx", "ScoreType", "Score")]
  names(eQTL.gtex1) <- c("SNP", "Gene", "Ref", "Tissue", "ScoreType", "Score")
 
  eQTL <- rbind(eQTL.uc, eQTL.gtex1)

  gwas.ld.eQTL <- merge(ld.snps, eQTL,  by.x = "SNP.ld",
                        by.y = "SNP", all.x = TRUE, all.y = FALSE)
  
  save(gwas.ld.eQTL, file = "./VarGeneMapping/gwas.ld.eQTL.RData", compress = TRUE) 
  
  invisible(gwas.ld.eQTL)  
  
}

# ld.fantom5: For each LD SNP, (1) indentify each fantom5 where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.fantom5 <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",  
                       fantom5.file="./FANTOM5/FANTOM5.RData",
                       mpos.rng.d = c("chr", "startpos", "endpos")){
  load(ld.file)
  load(fantom5.file)
  fantom5 <- FANTOM5
  names(fantom5)[names(fantom5)=="gene"] <- "Gene"
  fantom5$Gene <- as.character(fantom5$Gene)
  fantom5.d <- fantom5[,c("chr", "startpos", "endpos", "Gene", "tissue")]
  fantom5.d$chr <- gsub("chr", "", fantom5.d$chr)
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    # no promoter
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    
    # distal
    fantom5.d.i <- subset(fantom5.d, chr==i)
    
    ## Expand range as one row per position for fantom5.d.i
    fantom5.d.i$Length <- fantom5.d.i[, mpos.rng.d[3]] - fantom5.d.i[, mpos.rng.d[2]] + 1
    fantom5.d.i.exp <- data.frame(Position = as.vector(unlist(apply(fantom5.d.i[, mpos.rng.d[2:3]], 1,
                                                                    function(x) x[1]:x[2]))),
                                  Gene = rep(as.character(fantom5.d.i$Gene), fantom5.d.i$Length),
                                  Tissue = rep(fantom5.d.i$tissue, fantom5.d.i$Length))
    fantom5.d.i.exp <- fantom5.d.i.exp[!duplicated(fantom5.d.i.exp),]
    tmp <- merge(ld.snps.i, fantom5.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$Tissue)] <- "d"
    tmp <- tmp2[,c("CHR.ld","POS.ld", "SNP.ld", "Gene",  "Tissue", "type")] 
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  gwas.ld.fantom5 <- res
  save(gwas.ld.fantom5, file = "./VarGeneMapping/gwas.ld.fantom5.RData", compress = TRUE) 
  
  invisible(res)  
  
}


# ld.Posgene: For each LD SNP, (1) identify genes falling in the 5Kb windows (upstream 5kb, downstream 5kb) 
ld.Posgene <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                       Posgene.file="./gencode_v20_RefSeq37_1.RData",
                       win.size = 5000){
  cat("Start merging Posgene data...", "\n")
  load(ld.file)
  load(Posgene.file)
  names(gencode.refseq)[names(gencode.refseq)=="chrom"] <- "CHROM"
  names(gencode.refseq)[names(gencode.refseq)=="Gene"] <- "LOCUS"
  names(gencode.refseq)[names(gencode.refseq)=="length"] <- "Length"
  gencode.refseq$start <-  gencode.refseq$start - win.size 
  gencode.refseq$end <-  gencode.refseq$end + win.size  
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    
    # distal
    gencode.refseq.i <- subset(gencode.refseq, CHROM==i)
    
    ## Expand range as one row per position for gencode.refseq.i
    gencode.refseq.i$Length <- gencode.refseq.i[, "end"] - gencode.refseq.i[, "start"] + 1
    gencode.refseq.i.exp <- data.frame(Position = as.vector(unlist(apply(gencode.refseq.i[, c("start","end")], 1,
                                                                         function(x) x[1]:x[2]))),
                                       Gene = rep(as.character(gencode.refseq.i$LOCUS), gencode.refseq.i$Length))
    gencode.refseq.i.exp <- gencode.refseq.i.exp[!duplicated(gencode.refseq.i.exp),]
    tmp <- merge(ld.snps.i, gencode.refseq.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp <- tmp2[,c("CHR.ld","POS.ld", "SNP.ld", "Gene")] 
    
    cpos[[i]] <- tmp
  }
  gwas.ld.posgene <- do.call("rbind", cpos)
  save(gwas.ld.posgene, file = "./VarGeneMapping/gwas.ld.posgene.RData", compress = TRUE) 
  
  invisible(gwas.ld.posgene)  
  
}

# ld.vep: Merge gwas ld SNPs and VEP information. 
ld.vep <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                   vep.path="./1KG/Annotations"){
  cat("Start merging VEP annotation information...", "\n")
  load(ld.file)
  
  vep.files <- list.files(vep.path, pattern="*.txt.gz$") 
  n <- length(vep.files)
  chrs <- unlist(lapply( strsplit(vep.files,"_"),function(x)x[[3]]))
  chrs <- unlist(lapply( strsplit(chrs,"\\."),function(x)x[[1]]))
  chrs <- gsub("chr", "", chrs)
  m.no <- match(c(as.character(1:22), "X"), chrs)  
  
  gwas.ld.vep <- list() 
  k <- 0
  for (i in m.no){
    k <- k + 1
    cat(k, "\n")
    chr <- chrs[i]
    ld.snps.i <- subset(ld.snps, CHR.ld==chr)
    
    file <- vep.files[i]
    file <- paste(vep.path, file, sep="/")
    vep.data <- read.delim(file, skip=17,
                           header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
    names(vep.data)[names(vep.data)=="X.Uploaded_variation"] <- "Uploaded_variation"
    vep.data1 <- subset(vep.data, Uploaded_variation %in% ld.snps.i$SNP.ld | Existing_variation %in% ld.snps.i$SNP.ld)
    vep.data1 <- vep.data1[,c("Uploaded_variation", "SYMBOL",  "CANONICAL", "Consequence", "Amino_acids", "Condel", "Existing_variation")]
    
    ld.vep.d1 <- merge(ld.snps.i, vep.data1[,-ncol(vep.data1)], by.x = "SNP.ld",
                       by.y = "Uploaded_variation", all.x = TRUE, all.y = FALSE)
    ld.vep.d2 <- merge(ld.snps.i, vep.data1[,-1], by.x = "SNP.ld",
                       by.y = "Existing_variation", all.x = TRUE, all.y = FALSE)
    ld.vep.d <- merge(ld.vep.d1, ld.vep.d2, 
                      by = names(ld.vep.d1), 
                      all.x = TRUE, all.y = TRUE)
    
    gwas.ld.vep[[k]] <- ld.vep.d
  }
  gwas.ld.vep <- do.call("rbind", gwas.ld.vep)  
  gwas.ld.vep <- gwas.ld.vep[!duplicated(gwas.ld.vep),]
  save(gwas.ld.vep, file = "./VarGeneMapping/gwas.ld.vep.full.RData", compress = TRUE) 
  #cons <- sort(unique(gwas.ld.vep$Consequence))
  #writeLines(cons, "./VarGeneMapping/VEP_Consequence.txt", sep = "\n", useBytes = FALSE)
  
  invisible(gwas.ld.vep)  
  
}

# ld.vep.simplify: Simplify the full version of merged gwas ld SNPs and VEP data. 
ld.vep.simplify<- function(vep.file="./VarGeneMapping/gwas.ld.vep.full.RData", 
                           severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt"){
  cat("Start simplifying the merged VEP annotation data...", "\n")
  load(vep.file)
  severity <- read.delim(severity.file,
                         header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  
  snp.count <- count(gwas.ld.vep, c("SNP.ld", "SYMBOL"))
  snp.count1 <- subset(snp.count, freq == 1)
  snp.count2 <- subset(snp.count, freq > 1)
  snp.count1$id <- paste(snp.count1$SNP.ld, snp.count1$SYMBOL, sep="_")
  snp.count2$id <- paste(snp.count2$SNP.ld, snp.count2$SYMBOL, sep="_")
  
  gwas.ld.vep$id <- paste(gwas.ld.vep$SNP.ld, gwas.ld.vep$SYMBOL, sep="_")
  gwas.ld.vep1 <- subset(gwas.ld.vep, id %in% snp.count1$id)
  gwas.ld.vep2 <- subset(gwas.ld.vep, id %in% snp.count2$id)
  
  gwas.ld.vep3 <- subset(gwas.ld.vep2, CANONICAL %in% "YES")
  gwas.ld.vep4 <- subset(gwas.ld.vep2, !(CANONICAL %in% "YES"))
  gwas.ld.vep5 <- subset(gwas.ld.vep4, !(id %in%  gwas.ld.vep3$id))
  gwas.ld.vep5$severity <- severity$Severity[match(gwas.ld.vep5$Consequence, severity$Consequences)]
  gwas.ld.vep5 <- gwas.ld.vep5[order(gwas.ld.vep5$id, gwas.ld.vep5$severity),]
  gwas.ld.vep6 <- gwas.ld.vep5[match(unique(gwas.ld.vep5$id), gwas.ld.vep5$id),]
  
  gwas.ld.vep1$id <- NULL
  gwas.ld.vep3$id <- NULL
  gwas.ld.vep6$severity <- NULL
  gwas.ld.vep6$id <- NULL
  
  gwas.ld.vep.simplify <- rbind(gwas.ld.vep1, gwas.ld.vep3, gwas.ld.vep6)
  # sum(duplicated(gwas.ld.vep.simplify))  # 0
  save(gwas.ld.vep.simplify, file = "./VarGeneMapping/gwas.ld.vep.simplify.RData", compress = TRUE) 
  
  invisible(gwas.ld.vep.simplify)  
}


#merge.ld.simp: Merge based on SNP.ld and gene from each data
# % of unique (ld) GWAS SNPs with at least one genes
merge.ld.simp <- function(dhscor.file ="./VarGeneMapping/gwas.ld.rdb.dhscor.RData",
                     chiapet.file = "./VarGeneMapping/gwas.ld.chiapet.RData",
                     fantom5.file = "./VarGeneMapping/gwas.ld.fantom5.RData",
                     eQTL.file = "./VarGeneMapping/gwas.ld.eQTL.RData",
                     posgene.file = "./VarGeneMapping/gwas.ld.posgene.RData",
                     vep.file.simplified = "./VarGeneMapping/gwas.ld.vep.simplify.RData"){
  cat("Start merging the dhscor, chip-pet, fantom5 and vep data ...", "\n")
  load(dhscor.file)
  load(chiapet.file)
  load(fantom5.file)
  load(eQTL.file)
  load(posgene.file)
  load(vep.file.simplified)
  
  gwas.ld.rdb.dhscor <- rename(gwas.ld.rdb.dhscor, c(Cor = 'DHS.Cor',
                                                     type = "DHS.Type"))
  gwas.ld.chiapet <- rename(gwas.ld.chiapet, c(CellType = 'CHIAPET.CellType',
                                               type = "CHIAPET.Type"))
  gwas.ld.fantom5 <- rename(gwas.ld.fantom5, c(Tissue = 'FANTOM5.Tissue',
                                               type = "FANTOM5.Type"))
  gwas.ld.eQTL <- rename(gwas.ld.eQTL, c(Ref = 'eQTL.Ref',
                                         Score = 'eQTL.Score',
                                         ScoreType = "eQTL.ScoreType",
                                         Tissue = "eQTL.Tissue"))

  #gwas.ld.posgene <- rename(gwas.ld.posgene, c(Gene = 'Gene')) 
  gwas.ld.vep.simplify <- rename(gwas.ld.vep.simplify, c(SYMBOL = 'Gene', 
                                                         CANONICAL = 'VEP.Canonical',
                                                         Consequence = "VEP.Consequence",
                                                         Amino_acids = "VEP.AA",
                                                         Condel = "VEP.Condel"))
  
  gwas.ld.eQTL <- gwas.ld.eQTL[!duplicated(gwas.ld.eQTL),]
  # Remove the redundant rows (NA for Gene,  DHS.Cor, DHS.Type columns and their SNP.ld already presents) in the gwas.ld.rdb.dhscor data
  # gwas.ld.rdb.dhscor$id <- paste(gwas.ld.rdb.dhscor$Gene, gwas.ld.rdb.dhscor$DHS.Cor, gwas.ld.rdb.dhscor$DHS.Type, sep="_")
  # gwas.ld.rdb.dhscor1 <- subset(gwas.ld.rdb.dhscor, id!="NA_NA_NA") 
  # gwas.ld.rdb.dhscor2 <- subset(gwas.ld.rdb.dhscor, id=="NA_NA_NA") 
  # gwas.ld.rdb.dhscor3 <- subset(gwas.ld.rdb.dhscor2, !(SNP.ld %in% gwas.ld.rdb.dhscor1$SNP.ld))
  # gwas.ld.rdb.dhscor3$id <- NULL
  # gwas.ld.rdb.dhscor1$id <- NULL
  # gwas.ld.rdb.dhscor <- rbind(gwas.ld.rdb.dhscor1, gwas.ld.rdb.dhscor3)
  # gwas.ld.rdb.dhscor <- gwas.ld.rdb.dhscor[order(gwas.ld.rdb.dhscor$CHR.ld, gwas.ld.rdb.dhscor$POS.ld, gwas.ld.rdb.dhscor$SNP.ld),]
  #save(gwas.ld.rdb.dhscor, file = "./VarGeneMapping/GWAS_LD_SNPs_rdb_dhscor_updated.RData", compress = TRUE) 
  # Remove the redundant rows (NA for Gene, CHIAPET.CellType, CHIAPET.Type columns and their SNP.ld already presents) in the gwas.ld.chiapet data
  # gwas.ld.chiapet$id <- paste(gwas.ld.chiapet$Gene, gwas.ld.chiapet$CHIAPET.CellType, gwas.ld.chiapet$CHIAPET.Type, sep="_")
  # gwas.ld.chiapet1 <- subset(gwas.ld.chiapet, id!="NA_NA_NA") 
  # gwas.ld.chiapet2 <- subset(gwas.ld.chiapet, id=="NA_NA_NA") 
  # gwas.ld.chiapet3 <- subset(gwas.ld.chiapet2, !(SNP.ld %in% gwas.ld.chiapet1$SNP.ld))
  # gwas.ld.chiapet3$id <- NULL
  # gwas.ld.chiapet1$id <- NULL
  # gwas.ld.chiapet <- rbind(gwas.ld.chiapet1, gwas.ld.chiapet3)
  # gwas.ld.chiapet <- gwas.ld.chiapet[order(gwas.ld.chiapet$CHR.ld, gwas.ld.chiapet$POS.ld, gwas.ld.chiapet$SNP.ld),]
  #save(gwas.ld.chiapet, file = "./VarGeneMapping/GWAS_LD_SNPs_chiapet_updated.RData", compress = TRUE) 
  
  gwas.ld.rdb.dhscor1 <- subset(gwas.ld.rdb.dhscor, !is.na(Gene))
  gwas.ld.chiapet1 <- subset(gwas.ld.chiapet, !is.na(Gene))
  gwas.ld.fantom51 <- subset(gwas.ld.fantom5, !is.na(Gene))
  gwas.ld.eQTL1 <- subset(gwas.ld.eQTL, !is.na(Gene))
  gwas.ld.posgene1 <- subset(gwas.ld.posgene, !is.na(Gene))
  gwas.ld.vep.simplify1 <- subset(gwas.ld.vep.simplify, !is.na(Gene) & gwas.ld.vep.simplify$Gene!="-")
  
  gwas.ld.d <- merge(gwas.ld.rdb.dhscor1, gwas.ld.chiapet1[,-match(c("CHR.ld", "POS.ld"), names(gwas.ld.chiapet1))], 
                     by = c("SNP.ld", "Gene"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.fantom51[,-match(c("CHR.ld", "POS.ld"), names(gwas.ld.fantom51))], 
                     by = c("SNP.ld", "Gene"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.eQTL1[,-match(c("CHR.ld", "POS.ld"), names(gwas.ld.eQTL1))], 
                     by = c("SNP.ld", "Gene"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.posgene1[,-match(c("CHR.ld", "POS.ld"), names(gwas.ld.posgene1))], 
                     by = c("SNP.ld", "Gene"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.vep.simplify1[,-match(c("CHR.ld", "POS.ld"), names(gwas.ld.vep.simplify1))], 
                     by = c("SNP.ld", "Gene"), all = TRUE)
  
  # Add the CHR.ld, POS.ld back
  gwas.ld.d$CHR.ld <- gwas.ld.rdb.dhscor$CHR.ld[match(gwas.ld.d$SNP.ld, gwas.ld.rdb.dhscor$SNP.ld)]
  gwas.ld.d$POS.ld <- gwas.ld.rdb.dhscor$POS.ld[match(gwas.ld.d$SNP.ld, gwas.ld.rdb.dhscor$SNP.ld)]
  
  var2gene.simplified <- gwas.ld.d
  if(sum(duplicated(var2gene.simplified))>0){
    var2gene.simplified <- var2gene.simplified[!duplicated(var2gene.simplified),]
  }
  save(var2gene.simplified, file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData", compress = TRUE)
  invisible(var2gene.simplified)
}

# Subset of the var2gene data by using the union of RefSeq and Gencode gene names
# Also Update the var2gene mapping data by removing any rows with conflicting chr.ld with those from the chr information from gencode data
# 
var2gene.filter <- function(var2gene.file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData",
                            gencode.file = "./gencode_v20.RData",
                            RefSeq.file = "./RefSeq37_1.RData"){
  
  cat("Start to filter var2gene data ...", "\n")
  load(var2gene.file)
  #load(var2gene.file1)
  load(gencode.file)
  load(RefSeq.file)
  
  # Filter down the gene names not in the union of RefSeq and Gencode gene names
  # Combine the gencode and gene.data 
  # First filtering out any genes from RefSeq with more than one chromosome (they are transfer RNAs)
  gene.data1 <- gene.data[,c("LOCUS", "CHROM")]
  num <- count(gene.data1, c("LOCUS"))
  num1 <- subset(num, freq>1)
  gene.data1 <- subset(gene.data1, !(LOCUS %in% num1$LOCUS))
  gene.data1$source <- "refseq"
  
  gencode1 <- gencode[,c("Gene", "chrom")]
  gencode1$source <- "gencode"
  names(gene.data1) <- names(gencode1)
  
  gencode1 <- rbind(gencode1, gene.data1)
  gencode1 <- gencode1[!duplicated(gencode1[,-ncol(gencode1)]),]
  # Find the genes with more than 1 chrs
  num <- count(gencode1, c("Gene"))
  num1 <- subset(num, freq>1)
  # For genes in num1 (with different chr information in gencode and refseq), stick to the data from gencode by removing the data from RefSeq
  gencode1 <- subset(gencode1, !((Gene %in% num1$Gene) & (source %in% "refseq") ))
  
  genes <- unique(gencode1$Gene)
  var2gene.vepSimp <- subset(var2gene.simplified, Gene %in% genes)
  #var2gene.vepFull <- subset(var2gene.full, Gene %in% genes)
  
  # Also Update the var2gene mapping data by removing any rows with conflicting chr.ld with those from the chr information from gencode data
  var2gene.vepSimp$chr.gencode <- gencode1$chrom[match(var2gene.vepSimp$Gene, gencode1$Gene)]
  #var2gene.vepFull$chr.gencode <- gencode1$chrom[match(var2gene.vepFull$Gene, gencode1$Gene)]
  var2gene.vepSimp <- subset(var2gene.vepSimp, CHR.ld == chr.gencode)
  #var2gene.vepFull <- subset(var2gene.vepFull, CHR.ld == chr.gencode)  
  
  save(var2gene.vepSimp, file = "./VarGeneMapping/Var2Gene_vepSimp.RData", compress = TRUE)
  #save(var2gene.vepFull, file = "./VarGeneMapping/Var2Gene_vepFull.RData", compress = TRUE)
  invisible(var2gene.vepSimp)
}

# gwas.mesh: add mesh terms to the GWAS data. 
# Run it on Windows only
gwas.mesh <- function(gwas.file="./stopgap_3sources_dbSNP141.RData",
                      disease2mesh.file = "./disease.msh.final.txt"){
  
  cat("Start adding mesh terms to the GWAS data ...", "\n")
  load(gwas.file)
  dis2gene <- read.delim(disease2mesh.file,
                         header = TRUE, comment.char = "", na.strings = c("",".",NA), sep = "\t", as.is = TRUE) 
  
  gwas.data$disease <- tolower(trimWhiteSpace(gwas.data$disease))  
  gwas.data <- gwas.filter(gwas.data)
  dis2gene$disease <- tolower(dis2gene$disease)
  dis2gene$disease <- trimWhiteSpace(dis2gene$disease)
  
  #  all(unique(gwas.data$disease) %in% dis2gene$disease)  # FALSE, two disease names look so weird, remove them. 
  # unique(gwas.data$disease)[!(unique(gwas.data$disease) %in% unique(dis2gene$disease))]
  gwas.data <- subset(gwas.data, disease %in%  unique(dis2gene$disease))
  # on Unix, the tolower() will have error: Error in tolower(gwas.data$disease) : invalid multibyte string 3997
  # Run on Windows
  
  dis2gene$MSH <- tolower(dis2gene$MSH)
  dis2gene$MSH <- trimWhiteSpace(dis2gene$MSH)
  m.no <- match(gwas.data$disease, dis2gene$disease)
  gwas.data$msh <- dis2gene$MSH[m.no]
  gwas.data$msh.tree <- dis2gene$MSH_tree[m.no]
  gwas.data$msh.cat <- dis2gene$Category[m.no]
  keep.names <- c("snp_id","rsid","PUBMEDID", "pvalue","disease","msh", "msh.tree", "msh.cat")
  rest.names <- names(gwas.data)[!(names(gwas.data) %in% keep.names)]
  gwas.data <- gwas.data[,c(keep.names, rest.names)]
  save(gwas.data, file = "./stopgap_3sources_dbSNP141_msh.RData", compress = TRUE)  
  
  invisible(gwas.data)
}


# Merge GWAS data, LD r2 data and var2gene_vepSimp data together
gwas.ld.var2gene1 <- function(gwas.file="stopgap_3sources_dbSNP141_msh.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData",
                              var2gene.file = "./VarGeneMapping/Var2Gene_vepSimp.RData"){
  
  cat("Start merging GWAS data, LD r2 data and var2gene data together ...", "\n")
  load(gwas.file)
  load(ld.file)
  load(var2gene.file)
  var2gene.simplified <- var2gene.vepSimp
  var2gene.simplified$chr.gencode <- NULL
  
  gwas.r2 <- merge(gwas.data, ld.snps.r2[,c("SNP.gwas", "SNP.ld", "r2", "AF", "ASN_AF", "AMR_AF", "AFR_AF", "EUR_AF")], 
                   by.x = "rsid", by.y="SNP.gwas", all = TRUE)
  gwas.r2.var2gene.full <- merge(gwas.r2, var2gene.simplified, 
                                 by.x = "SNP.ld", by.y="SNP.ld", all = TRUE)
  gwas.r2.var2gene.withGene <- subset(gwas.r2.var2gene.full, !is.na(Gene))
  new.names <- c("snp.ld", "snp.gwas", "snp.gwas.orig", "pubmedid", "pvalue",
                 "disease","msh","msh.tree","msh.cat", "init.samp", "source", "rep.samp", "nhlbikey", 
                 "loc.paper", "datepub", "gwas.ancestry", "sampsize.dis.rep",
                 "sampsize.dis", "sampsize.rep", "snp.gwasdb", "source.gwasdb",
                 "hpo.id", "hpo.term", "do.id", "do.term", "r2", "af.1kg", "asn.af.1kg", 
                 "amr.af.1kg", "afr.af.1kg", "eur.af.1kg","gene", "chr.ld",
                 "pos.ld", "cat.rdb", "dhs.cor", "dhs.type", "chiapet.cell", "chiapet.type",
                 "fantom5.tissue", "fantom5.type", "eqtl.ref","eqtl.tissue","eqtl.scoretype",
                 "eqtl.score", "vep.canonical", "vep.conseq",
                 "vep.aa", "vep.condel")
  names(gwas.r2.var2gene.full) <- new.names
  names(gwas.r2.var2gene.withGene) <- new.names
  
  # Remove  two papers with lots of metabolite associations: PUBMEDID = 23281178 and 22286219
  # Remove papers with large numbers of associations (>500) and identified these papers to be excluded.
  #   PUBMEDID  disease  Exclude
  #  17903293  Biomarkers	TRUE
  #  19043545	Metabolite levels	TRUE
  #  22286219	Metabolite levels	TRUE
  #  23281178	Metabolite levels	TRUE
  #  17554300	Multiple complex diseases	TRUE
  #  19862010	Other erythrocyte phenotypes	TRUE
  gwas.r2.var2gene.full<- subset(gwas.r2.var2gene.full, !(pubmedid %in% c("23281178","22286219", "17903293", "19043545", "17554300", "19862010")))
  gwas.r2.var2gene.withGene<- subset(gwas.r2.var2gene.withGene, !(pubmedid %in% c("23281178","22286219", "17903293", "19043545", "17554300", "19862010")))                                
  
  system("mkdir GWAS_LD_var2gene")
  save(gwas.r2.var2gene.full, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_full_r2_0.5.RData", compress = TRUE)
  save(gwas.r2.var2gene.withGene, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5.RData", compress = TRUE)
  
  gwas.r2.var2gene.full1 <- subset(gwas.r2.var2gene.full, r2>=0.7)
  save(gwas.r2.var2gene.full1, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_full_r2_0.7.RData", compress = TRUE)
  gwas.r2.var2gene.withGene1 <- subset(gwas.r2.var2gene.withGene, r2>=0.7)
  save(gwas.r2.var2gene.withGene1, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7.RData", compress = TRUE)
  
  invisible(gwas.r2.var2gene.full)
  
}

# find MHC genes:
find.mhc.genes <- function(gencode.file = "./gencode_v20_RefSeq37_1.RData")
{
  load(gencode.file)
  x <- subset(gencode.refseq, chrom %in% "6" & start >= 26016335 & end <= 33548071)
  mhc <- unique(x$Gene)
  
  invisible(mhc)
}

# remove MHC genes from the merged r2=0.7 data:
rm.mhc.genes <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7.RData",
                         lRmMHC=TRUE){
  
  load(data.file)
  mhc.genes <- find.mhc.genes(gencode.file = "./gencode_v20_RefSeq37_1.RData")
  write.table(mhc.genes, "./GWAS_LD_var2gene/MHC_genes.txt", sep = "\t", row.names = F)
  if (lRmMHC){
    gwas.r2.var2gene.withGene.rmMHC <- subset(gwas.r2.var2gene.withGene1, !(gene %in% mhc.genes))  
  } else {
    gwas.r2.var2gene.withGene.rmMHC <- gwas.r2.var2gene.withGene1
  }
  
  save(gwas.r2.var2gene.withGene.rmMHC, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC.RData", compress = TRUE)
  
  invisible(gwas.r2.var2gene.withGene.rmMHC)
}


# Add gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC.RData data
gene.score <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC.RData",
                       severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt"){
  
  cat("Start adding gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7.RData data ...", "\n")
  load(data.file)
  severity <- read.delim(severity.file,
                         header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  gwas.r2.var2gene.withGene.rmMHC$severity <- severity$Severity[match(gwas.r2.var2gene.withGene.rmMHC$vep.conseq, severity$Consequences)]
  
  # 1: Source
  s.lookup <- data.frame(Source=c("nhgri", "grasp", "gwasdb"), score=c(3,2,2))
  gwas.r2.var2gene.withGene.rmMHC$score1 <- s.lookup$score[match(gwas.r2.var2gene.withGene.rmMHC$source, s.lookup$Source)]
  
  # 2: r2
  gwas.r2.var2gene.withGene.rmMHC$r2.bin <- as.character(cut(gwas.r2.var2gene.withGene.rmMHC$r2, breaks=c(0.7,0.8,0.9,1), right=FALSE, include.highest=TRUE))
  gwas.r2.var2gene.withGene.rmMHC$r2.bin[which((gwas.r2.var2gene.withGene.rmMHC$r2) <= 1 & (gwas.r2.var2gene.withGene.rmMHC$r2) >= 0.9)] <- "[0.9,1]"
  r2.lookup <- data.frame(r2.bin=c("[0.7,0.8)", "[0.8,0.9)", "[0.9,1]"), score=c(0,1,2))
  gwas.r2.var2gene.withGene.rmMHC$score2 <- r2.lookup$score[match(gwas.r2.var2gene.withGene.rmMHC$r2.bin, r2.lookup$r2.bin)]
  
  # 3: cat.rdb
  cat.lookup <- data.frame(cat=c("1a","1b","1c","1d","1e","1f","2a","2b","2c", "3a","3b", "4"), score=c(rep(2, 9), rep(1, 2), 0))
  gwas.r2.var2gene.withGene.rmMHC$score3 <- cat.lookup$score[match(gwas.r2.var2gene.withGene.rmMHC$cat.rdb, cat.lookup$cat)]
  
  # 3.1: dhs.cor  # Add this one as well
  gwas.r2.var2gene.withGene.rmMHC$score3.1 <- NA
  gwas.r2.var2gene.withGene.rmMHC$score3.1[!is.na(gwas.r2.var2gene.withGene.rmMHC$dhs.type)] <- 2
  
  # 4: chiapet
  gwas.r2.var2gene.withGene.rmMHC$score4 <- NA
  gwas.r2.var2gene.withGene.rmMHC$score4[!is.na(gwas.r2.var2gene.withGene.rmMHC$chiapet.cell)] <- 1
  # if presenting in multiple cell lines, further add 1 
  gwas.r2.var2gene.withGene.rmMHC$id <- apply(gwas.r2.var2gene.withGene.rmMHC[,c("snp.ld", "snp.gwas", "gene")], 1, paste, collapse="_")
  gwas.r2.var2gene.withGene.rmMHC1 <- subset(gwas.r2.var2gene.withGene.rmMHC, !is.na(chiapet.cell))

  id.count <- count(gwas.r2.var2gene.withGene.rmMHC1, "id")
  id.count1 <- subset(id.count, freq>1) # including ids presenting in multiple cell lines
  gwas.r2.var2gene.withGene.rmMHC$score4[gwas.r2.var2gene.withGene.rmMHC$id %in% id.count1] <- 2
  
  # 5: fantom5
  gwas.r2.var2gene.withGene.rmMHC$score5.1 <- gwas.r2.var2gene.withGene.rmMHC$score5.2 <- NA
  gwas.r2.var2gene.withGene.rmMHC$score5.1[!is.na(gwas.r2.var2gene.withGene.rmMHC$fantom5.tissue)] <- 1
  gwas.r2.var2gene.withGene.rmMHC$score5.2[grep(", ", gwas.r2.var2gene.withGene.rmMHC$fantom5.tissue)] <- 1
  gwas.r2.var2gene.withGene.rmMHC$score5 <- apply(gwas.r2.var2gene.withGene.rmMHC[,c("score5.1", "score5.2")], 1, sum, na.rm=TRUE)
  gwas.r2.var2gene.withGene.rmMHC$score5.1 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score5.2 <- NULL
  
  # 6: eqtl
  gwas.r2.var2gene.withGene.rmMHC$score6 <- NA
  gwas.r2.var2gene.withGene.rmMHC$score6[!is.na(gwas.r2.var2gene.withGene.rmMHC$eqtl.ref)] <- 2
  # if presenting in multiple cell lines, further add 1 
  gwas.r2.var2gene.withGene.rmMHC1 <- subset(gwas.r2.var2gene.withGene.rmMHC, !is.na(eqtl.ref))
  id.count <- count(gwas.r2.var2gene.withGene.rmMHC1, "id")
  id.count1 <- subset(id.count, freq>1)  # including ids presenting in multiple studies
  gwas.r2.var2gene.withGene.rmMHC$score6[gwas.r2.var2gene.withGene.rmMHC$id %in% id.count1$id] <- 3
  rm(gwas.r2.var2gene.withGene.rmMHC1)
  
  # 7: vep
  gwas.r2.var2gene.withGene.rmMHC$score7 <- NA
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity <= 5] <- 4
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity <= 10 & gwas.r2.var2gene.withGene.rmMHC$severity >= 6] <- 3
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity <= 17 & gwas.r2.var2gene.withGene.rmMHC$severity >= 11] <- 2
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity <= 33 & gwas.r2.var2gene.withGene.rmMHC$severity >= 18] <- 1
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity == 34] <- 0
  # Severity 6~10: for missense SNPs: +3 if condel deleterious, else neutral +2
  is.neu <- rep(FALSE, nrow(gwas.r2.var2gene.withGene.rmMHC))
  is.neu[grep("neutral", gwas.r2.var2gene.withGene.rmMHC$vep.condel)] <- TRUE
  gwas.r2.var2gene.withGene.rmMHC$score7[gwas.r2.var2gene.withGene.rmMHC$severity <= 10 & 
                                           gwas.r2.var2gene.withGene.rmMHC$severity >= 6 & 
                                           gwas.r2.var2gene.withGene.rmMHC$vep.conseq =="missense_variant" & 
                                           is.neu] <- 2
  
  gwas.r2.var2gene.withGene.rmMHC$gene.score <- apply(gwas.r2.var2gene.withGene.rmMHC[,c("score1", "score2", "score3", "score3.1", "score4", "score5", "score6", "score7")], 1, sum, na.rm=TRUE)
  gwas.r2.var2gene.withGene.rmMHC$score1 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score2 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score3 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score3.1 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score4 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score5 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score6 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$score7 <- NULL
  gwas.r2.var2gene.withGene.rmMHC$r2.bin <- NULL
  gwas.r2.var2gene.withGene.rmMHC$id <- NULL
  save(gwas.r2.var2gene.withGene.rmMHC, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore.RData", compress = TRUE)
  
  invisible(gwas.r2.var2gene.withGene.rmMHC)
  
}

# For each unique GWAS SNP, rank its LD SNPs' genes based on their gene scores and r2 if gene scores tie
gene.rank <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore.RData"){
  
  cat("For each unique GWAS SNP, start ranking its LD SNPs' genes based on their gene scores...", "\n")
  load(data.file)
  data <- gwas.r2.var2gene.withGene.rmMHC[order(gwas.r2.var2gene.withGene.rmMHC$snp.gwas,
                                                -gwas.r2.var2gene.withGene.rmMHC$gene.score,
                                                -gwas.r2.var2gene.withGene.rmMHC$r2),]
  data$id <- paste(data$snp.gwas, data$gene, sep="_")
  #data.a <- subset(data, snp.gwas %in% c("rs13111494","rs884304", "rs4665058"))
  data.a <- data[, c("snp.gwas", "gene")]
  data.a1 <- data.a[!duplicated(data.a),]
  data.a1$id <- paste(data.a1$snp.gwas, data.a1$gene, sep="_")
  
  count1 <- count(data.a, c("snp.gwas", "gene"))  
  # Please note count() automatically sorts the snp.gwas-gene pair names, so need to put the sorted order back
  count2 <- count(count1[,-3], "snp.gwas")    # Need to remove the freq column (-3), otherwise the numbers are not what I want
  count1$id <- paste(count1$snp.gwas, count1$gene,sep="_")
  count1 <- count1[match(data.a1$id, count1$id),]
  count1$rank <- unlist(lapply(count2$freq, function(x) 1:x))
  #data.a$gene <- as.character(data.a$gene)
  #data$snp.gwas <- as.character(data$snp.gwas)
  
  data$gene.rank.random <- count1$rank[match(data$id, count1$id)]
  
  # Rank based on rank(x, ties.method = "max") and rank(x, ties.method = "min")
  data1 <- data[,c("snp.gwas", "gene", "id", "gene.score", "r2")]
  data2 <- data1[match(count1$id, data1$id),] 
  data2$genescorer2 <- -1 * data2$gene.score * data2$r2
  # rank 1 means the best gene score and r2
  gene.rank1 <- aggregate(genescorer2 ~ snp.gwas, data2, function(x){rank(x, ties.method = "max")})
  # Divide out ranks with more than one rank 
  rank.max <- as.numeric(unlist(gene.rank1$genescorer2))
  count1$rank.max <- rank.max  
  data$gene.rank.max <- count1$rank.max[match(data$id, count1$id)]
  
  gene.rank2 <- aggregate(genescorer2 ~ snp.gwas, data2, function(x){rank(x, ties.method = "min")})
  # Divide out ranks with more than one rank 
  rank.min <- as.numeric(unlist(gene.rank2$genescorer2))
  count1$rank.min <- rank.min  
  data$gene.rank.min <- count1$rank.min[match(data$id, count1$id)]  
  
  # Remove the gene.rank.random column
  data$gene.rank.random <- NULL
  data$id <- NULL
  stopgap.data <- data
  save(stopgap.data, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank.RData", compress = TRUE)
  
  # Select the best ld SNPs per gene and each original GWAS SNP based their gene score and r2
  # Order by the snp.gwas, rev(gene.score) and rev(r2), then pick the first row
  
  data.bestld <- stopgap.data
  data.bestld$id <- paste(data.bestld$snp.gwas, data.bestld$pubmedid, data.bestld$disease, data.bestld$gene, sep="_")
  data.bestld <- data.bestld[match(unique(data.bestld$id), data.bestld$id),]
  data.bestld$id <- NULL
  
  stopgap.data.bestld <- data.bestld
  save(stopgap.data.bestld, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD.RData", compress = TRUE)
  
  invisible(data)
}


#ld.data: Get the LD data as a starting point for SNP locus clustering:
ld.data <- function(data.file="./LocusClustering/gwas_r2.RData"){
  load(data.file)
  gwas.r2.var2gene.full.ld <- gwas.r2[,c("snp.ld","snp.gwas","pvalue","r2","chr.ld", "pos.ld")]
  gwas.r2.var2gene.full.ld <- gwas.r2.var2gene.full.ld[!duplicated(gwas.r2.var2gene.full.ld),]
  gwas.r2.var2gene.full.ld <- gwas.r2.var2gene.full.ld[order(gwas.r2.var2gene.full.ld$pvalue,
                                                             gwas.r2.var2gene.full.ld$snp.gwas,
                                                             gwas.r2.var2gene.full.ld$snp.ld),]
  ld1 <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% c("1", "2", "3"))
  ld2 <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% c("4", "5", "6"))
  ld3 <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% c("7", "8", "9","10", "11"))
  ld4 <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% c("12", "13", "14","15", "16", "17"))
  ld5 <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% c("18","19", "20", "21", "22", "X", "Y"))
  #save(gwas.r2.var2gene.full.ld, file = "./LocusClustering/ld.RData", compress = TRUE)
  save(ld1, file = "./LocusClustering/ld1.RData", compress = TRUE)
  save(ld2, file = "./LocusClustering/ld2.RData", compress = TRUE)
  save(ld3, file = "./LocusClustering/ld3.RData", compress = TRUE)
  save(ld4, file = "./LocusClustering/ld4.RData", compress = TRUE)
  save(ld5, file = "./LocusClustering/ld5.RData", compress = TRUE)
  
  invisible(gwas.r2.var2gene.full.ld)
}

# locus.cluster1: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm does not allow the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster1 <- function(ld.file="./LocusClustering/ld1.RData"){
  load(ld.file)
  
  gwas.r2.var2gene.full.r2 <- ld1
  #gwas.r2.var2gene.full.r2 <- subset(ld1, snp.gwas %in% unique(ld1$snp.gwas)[1:10])
  
  gwas.r2.var2gene.full.r2$snp.gwas <- as.character(gwas.r2.var2gene.full.r2$snp.gwas)
  gwas.r2.var2gene.full.r2$snp.ld <- as.character(gwas.r2.var2gene.full.r2$snp.ld)
  gwas.snps <- unique(gwas.r2.var2gene.full.r2$snp.gwas)
  ld.snps <- unique(gwas.r2.var2gene.full.r2$snp.ld)
  
  gwas.wait.snps <- gwas.snps
  ld.wait.snps <- ld.snps
  wait.data <- gwas.r2.var2gene.full.r2
  i <- 0
  clusters <- list()
  while(length(gwas.wait.snps) > 0){
    i <- i + 1
    cat(i, "\n")
    
    gwas.snps.i <- gwas.wait.snps[1]
    start.gwassnp <- gwas.snps.i
    cluster.i <- subset(wait.data, snp.gwas %in% gwas.snps.i)
    wait.data <- subset(wait.data, !(snp.gwas %in% gwas.snps.i))
    
    #gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% wait.data$snp.gwas])
    if(nrow(wait.data) > 0){
      ld.snps.i <- unique(cluster.i$snp.ld)
      cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      
      while(nrow(cluster.i.d)>0){
        # Update the wait.data and cluster.i and ld.snps.i here:
        cluster.i <- rbind(cluster.i, cluster.i.d)
        wait.data <- subset(wait.data, !((snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i)))
        ld.snps.i <- unique(cluster.i$snp.ld)
        
        cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      }
      
      gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% gwas.r2.var2gene.full.r2$snp.gwas])
      gwas.wait.snps <- gwas.wait.snps[!(gwas.wait.snps %in% gwas.snps.i)]
      
      ld.wait.snps <- ld.wait.snps[!(ld.wait.snps %in% ld.snps.i)]
      cat(length(gwas.wait.snps), "\n")
      
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      
      clusters[[i]] <- cluster.i
    } else {
      gwas.wait.snps <- NULL
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      clusters[[i]] <- cluster.i
    }
  }
  ld1.clusters <- do.call("rbind", clusters)
  save(ld1.clusters, file = "./LocusClustering/ld1_LocusClusters.RData", compress = TRUE)
  invisible(ld1.clusters)
}

# locus.cluster2: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm allows the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster2 <- function(ld.file="./LocusClustering/ld2.RData"){
  load(ld.file)
  
  gwas.r2.var2gene.full.r2 <- ld2
  #gwas.r2.var2gene.full.r2 <- subset(ld2, snp.gwas %in% unique(ld2$snp.gwas)[1:10])
  
  gwas.r2.var2gene.full.r2$snp.gwas <- as.character(gwas.r2.var2gene.full.r2$snp.gwas)
  gwas.r2.var2gene.full.r2$snp.ld <- as.character(gwas.r2.var2gene.full.r2$snp.ld)
  gwas.snps <- unique(gwas.r2.var2gene.full.r2$snp.gwas)
  ld.snps <- unique(gwas.r2.var2gene.full.r2$snp.ld)
  
  gwas.wait.snps <- gwas.snps
  ld.wait.snps <- ld.snps
  wait.data <- gwas.r2.var2gene.full.r2
  i <- 0
  clusters <- list()
  while(length(gwas.wait.snps) > 0){
    i <- i + 1
    cat(i, "\n")
    
    gwas.snps.i <- gwas.wait.snps[1]
    start.gwassnp <- gwas.snps.i
    cluster.i <- subset(wait.data, snp.gwas %in% gwas.snps.i)
    wait.data <- subset(wait.data, !(snp.gwas %in% gwas.snps.i))
    
    #gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% wait.data$snp.gwas])
    if(nrow(wait.data) > 0){
      ld.snps.i <- unique(cluster.i$snp.ld)
      cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      
      while(nrow(cluster.i.d)>0){
        # Update the wait.data and cluster.i and ld.snps.i here:
        cluster.i <- rbind(cluster.i, cluster.i.d)
        wait.data <- subset(wait.data, !((snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i)))
        ld.snps.i <- unique(cluster.i$snp.ld)
        
        cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      }
      
      gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% gwas.r2.var2gene.full.r2$snp.gwas])
      gwas.wait.snps <- gwas.wait.snps[!(gwas.wait.snps %in% gwas.snps.i)]
      
      ld.wait.snps <- ld.wait.snps[!(ld.wait.snps %in% ld.snps.i)]
      cat(length(gwas.wait.snps), "\n")
      
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      
      clusters[[i]] <- cluster.i
    } else {
      gwas.wait.snps <- NULL
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      clusters[[i]] <- cluster.i
    }
  }
  ld2.clusters <- do.call("rbind", clusters)
  save(ld2.clusters, file = "./LocusClustering/ld2_LocusClusters.RData", compress = TRUE)
  invisible(ld2.clusters)
}

# locus.cluster3: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm allows the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster3 <- function(ld.file="./LocusClustering/ld3.RData"){
  load(ld.file)
  
  gwas.r2.var2gene.full.r2 <- ld3
  #gwas.r2.var2gene.full.r2 <- subset(ld3, snp.gwas %in% unique(ld3$snp.gwas)[1:10])
  
  gwas.r2.var2gene.full.r2$snp.gwas <- as.character(gwas.r2.var2gene.full.r2$snp.gwas)
  gwas.r2.var2gene.full.r2$snp.ld <- as.character(gwas.r2.var2gene.full.r2$snp.ld)
  gwas.snps <- unique(gwas.r2.var2gene.full.r2$snp.gwas)
  ld.snps <- unique(gwas.r2.var2gene.full.r2$snp.ld)
  
  gwas.wait.snps <- gwas.snps
  ld.wait.snps <- ld.snps
  wait.data <- gwas.r2.var2gene.full.r2
  i <- 0
  clusters <- list()
  while(length(gwas.wait.snps) > 0){
    i <- i + 1
    cat(i, "\n")
    
    gwas.snps.i <- gwas.wait.snps[1]
    start.gwassnp <- gwas.snps.i
    cluster.i <- subset(wait.data, snp.gwas %in% gwas.snps.i)
    wait.data <- subset(wait.data, !(snp.gwas %in% gwas.snps.i))
    
    #gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% wait.data$snp.gwas])
    if(nrow(wait.data) > 0){
      ld.snps.i <- unique(cluster.i$snp.ld)
      cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      
      while(nrow(cluster.i.d)>0){
        # Update the wait.data and cluster.i and ld.snps.i here:
        cluster.i <- rbind(cluster.i, cluster.i.d)
        wait.data <- subset(wait.data, !((snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i)))
        ld.snps.i <- unique(cluster.i$snp.ld)
        
        cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      }
      
      gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% gwas.r2.var2gene.full.r2$snp.gwas])
      gwas.wait.snps <- gwas.wait.snps[!(gwas.wait.snps %in% gwas.snps.i)]
      
      ld.wait.snps <- ld.wait.snps[!(ld.wait.snps %in% ld.snps.i)]
      cat(length(gwas.wait.snps), "\n")
      
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      
      clusters[[i]] <- cluster.i
    } else {
      gwas.wait.snps <- NULL
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      clusters[[i]] <- cluster.i
    }
  }
  ld3.clusters <- do.call("rbind", clusters)
  save(ld3.clusters, file = "./LocusClustering/ld3_LocusClusters.RData", compress = TRUE)
  invisible(ld3.clusters)
}

# locus.cluster4: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm allows the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster4 <- function(ld.file="./LocusClustering/ld4.RData"){
  load(ld.file)
  
  gwas.r2.var2gene.full.r2 <- ld4
  #gwas.r2.var2gene.full.r2 <- subset(ld4, snp.gwas %in% unique(ld4$snp.gwas)[1:10])
  
  gwas.r2.var2gene.full.r2$snp.gwas <- as.character(gwas.r2.var2gene.full.r2$snp.gwas)
  gwas.r2.var2gene.full.r2$snp.ld <- as.character(gwas.r2.var2gene.full.r2$snp.ld)
  gwas.snps <- unique(gwas.r2.var2gene.full.r2$snp.gwas)
  ld.snps <- unique(gwas.r2.var2gene.full.r2$snp.ld)
  
  gwas.wait.snps <- gwas.snps
  ld.wait.snps <- ld.snps
  wait.data <- gwas.r2.var2gene.full.r2
  i <- 0
  clusters <- list()
  while(length(gwas.wait.snps) > 0){
    i <- i + 1
    cat(i, "\n")
    
    gwas.snps.i <- gwas.wait.snps[1]
    start.gwassnp <- gwas.snps.i
    cluster.i <- subset(wait.data, snp.gwas %in% gwas.snps.i)
    wait.data <- subset(wait.data, !(snp.gwas %in% gwas.snps.i))
    
    #gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% wait.data$snp.gwas])
    if(nrow(wait.data) > 0){
      ld.snps.i <- unique(cluster.i$snp.ld)
      cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      
      while(nrow(cluster.i.d)>0){
        # Update the wait.data and cluster.i and ld.snps.i here:
        cluster.i <- rbind(cluster.i, cluster.i.d)
        wait.data <- subset(wait.data, !((snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i)))
        ld.snps.i <- unique(cluster.i$snp.ld)
        
        cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      }
      
      gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% gwas.r2.var2gene.full.r2$snp.gwas])
      gwas.wait.snps <- gwas.wait.snps[!(gwas.wait.snps %in% gwas.snps.i)]
      
      ld.wait.snps <- ld.wait.snps[!(ld.wait.snps %in% ld.snps.i)]
      cat(length(gwas.wait.snps), "\n")
      
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      
      clusters[[i]] <- cluster.i
    } else {
      gwas.wait.snps <- NULL
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      clusters[[i]] <- cluster.i
    }
  }
  ld4.clusters <- do.call("rbind", clusters)
  save(ld4.clusters, file = "./LocusClustering/ld4_LocusClusters.RData", compress = TRUE)
  invisible(ld4.clusters)
}

# locus.cluster5: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm allows the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster5 <- function(ld.file="./LocusClustering/ld5.RData"){
  load(ld.file)
  
  gwas.r2.var2gene.full.r2 <- ld5
  #gwas.r2.var2gene.full.r2 <- subset(ld5, snp.gwas %in% unique(ld5$snp.gwas)[1:10])
  
  gwas.r2.var2gene.full.r2$snp.gwas <- as.character(gwas.r2.var2gene.full.r2$snp.gwas)
  gwas.r2.var2gene.full.r2$snp.ld <- as.character(gwas.r2.var2gene.full.r2$snp.ld)
  gwas.snps <- unique(gwas.r2.var2gene.full.r2$snp.gwas)
  ld.snps <- unique(gwas.r2.var2gene.full.r2$snp.ld)
  
  gwas.wait.snps <- gwas.snps
  ld.wait.snps <- ld.snps
  wait.data <- gwas.r2.var2gene.full.r2
  i <- 0
  clusters <- list()
  while(length(gwas.wait.snps) > 0){
    i <- i + 1
    cat(i, "\n")
    
    gwas.snps.i <- gwas.wait.snps[1]
    start.gwassnp <- gwas.snps.i
    cluster.i <- subset(wait.data, snp.gwas %in% gwas.snps.i)
    wait.data <- subset(wait.data, !(snp.gwas %in% gwas.snps.i))
    
    #gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% wait.data$snp.gwas])
    if(nrow(wait.data) > 0){
      ld.snps.i <- unique(cluster.i$snp.ld)
      cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      
      while(nrow(cluster.i.d)>0){
        # Update the wait.data and cluster.i and ld.snps.i here:
        cluster.i <- rbind(cluster.i, cluster.i.d)
        wait.data <- subset(wait.data, !((snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i)))
        ld.snps.i <- unique(cluster.i$snp.ld)
        
        cluster.i.d <- subset(wait.data, (snp.gwas %in% ld.snps.i) | (snp.ld %in% ld.snps.i))
      }
      
      gwas.snps.i <- unique(cluster.i$snp.ld[cluster.i$snp.ld %in% gwas.r2.var2gene.full.r2$snp.gwas])
      gwas.wait.snps <- gwas.wait.snps[!(gwas.wait.snps %in% gwas.snps.i)]
      
      ld.wait.snps <- ld.wait.snps[!(ld.wait.snps %in% ld.snps.i)]
      cat(length(gwas.wait.snps), "\n")
      
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      
      clusters[[i]] <- cluster.i
    } else {
      gwas.wait.snps <- NULL
      cluster.i$start.gwassnp <- start.gwassnp
      cluster.i$cluster <- i
      clusters[[i]] <- cluster.i
    }
  }
  ld5.clusters <- do.call("rbind", clusters)
  save(ld5.clusters, file = "./LocusClustering/ld5_LocusClusters.RData", compress = TRUE)
  invisible(ld5.clusters)
}

# cluster.merge: merge five clusters from different chrs together
cluster.merge <- function(c.file1="./LocusClustering/ld1_LocusClusters.RData",
                          c.file2="./LocusClustering/ld2_LocusClusters.RData",
                          c.file3="./LocusClustering/ld3_LocusClusters.RData",
                          c.file4="./LocusClustering/ld4_LocusClusters.RData",
                          c.file5="./LocusClustering/ld5_LocusClusters.RData"){
  load(c.file1)
  load(c.file2)
  load(c.file3)
  load(c.file4)
  load(c.file5)
  nc1 <- max(ld1.clusters$cluster)
  nc2 <- max(ld2.clusters$cluster)
  nc3 <- max(ld3.clusters$cluster)
  nc4 <- max(ld4.clusters$cluster)
  nc5 <- max(ld5.clusters$cluster)
  ld2.clusters$cluster <- ld2.clusters$cluster + nc1
  ld3.clusters$cluster <- ld3.clusters$cluster + nc1 + nc2
  ld4.clusters$cluster <- ld4.clusters$cluster + nc1 + nc2 + nc3
  ld5.clusters$cluster <- ld5.clusters$cluster + nc1 + nc2 + nc3 + nc4
  stopgap.snp.clusters <- rbind(ld1.clusters, ld2.clusters, ld3.clusters, ld4.clusters, ld5.clusters)
  save(stopgap.snp.clusters, file = "./LocusClustering/STOPGAP_SNP_Clusters.RData", compress = TRUE)
  
  invisible(stopgap.snp.clusters)
}

# Add cluser number data to the best LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank.RData",
                        cluster.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData"){
  
  cat("Start merging cluster data and the merged data before bestLD data together ...", "\n")
  load(data.file)
  load(cluster.file)
  
  #   stopgap.data$id <- paste(stopgap.data$snp.ld, stopgap.data$snp.gwas, sep="_")
  #   stopgap.snp.clusters$id <- paste(stopgap.snp.clusters$snp.ld, stopgap.snp.clusters$snp.gwas, sep="_")
  #   all(unique(stopgap.data$id) %in% unique(stopgap.snp.clusters$id))
  #   stopgap.data$id <- NULL
  #   stopgap.snp.clusters$id <- NULL
  
  stopgap.snp.clusters1 <- stopgap.snp.clusters
  names(stopgap.snp.clusters1)[names(stopgap.snp.clusters1)=="start.gwassnp"] <- "clust.init.gwassnp"
  
  stopgap.data <-   merge(stopgap.data, stopgap.snp.clusters1,
                          by = c("snp.ld","snp.gwas","chr.ld", "pos.ld", "pvalue","r2"), all.x = TRUE, all.y=FALSE)
  #sum(is.na(stopgap.data$cluster))
  
  # Rank the p-values of snp.gwas within each cluster in each pubmedid/disease combination
  stopgap.data <- stopgap.data[order(stopgap.data$pubmedid,
                                     stopgap.data$disease,
                                     stopgap.data$cluster,
                                     stopgap.data$pvalue),]
  stopgap.data$id <- paste(stopgap.data$pubmedid, stopgap.data$disease, stopgap.data$cluster, sep="_")
  stopgap.data$id1 <- paste(stopgap.data$id, stopgap.data$pvalue, sep="_")
  
  # Add the p-value rank
  tmp <- stopgap.data[,c("id", "id1")]
  tmp1 <- tmp[match(unique(tmp$id1),tmp$id1),]
  tmp2 <- tmp1[match(unique(tmp1$id),tmp1$id),]
  library(plyr)
  id.count <- count(tmp1, "id")
  # Put back the order. 
  id.count <- id.count[match(tmp2$id, id.count$id),]
  tmp1$rank <- unlist(lapply(id.count$freq, function(x) 1:x))
  stopgap.data$gwassnp.p.rank <- tmp1$rank[match(stopgap.data$id1, tmp1$id1)]
  stopgap.data$id <- NULL
  stopgap.data$id1 <- NULL
  rm(tmp)
  rm(tmp1)
  
  save(stopgap.data, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_cluster.RData", compress = TRUE)
  
}


# Add cluser number data to the best LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster.bestld <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD.RData",
                        cluster.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData"){
  
  cat("Start merging cluster data and bestLD data together ...", "\n")
  load(data.file)
  load(cluster.file)
  
  #   stopgap.data.bestld$id <- paste(stopgap.data.bestld$snp.ld, stopgap.data.bestld$snp.gwas, sep="_")
  #   stopgap.snp.clusters$id <- paste(stopgap.snp.clusters$snp.ld, stopgap.snp.clusters$snp.gwas, sep="_")
  #   all(unique(stopgap.data.bestld$id) %in% unique(stopgap.snp.clusters$id))
  #   stopgap.data.bestld$id <- NULL
  #   stopgap.snp.clusters$id <- NULL
  
  stopgap.snp.clusters1 <- stopgap.snp.clusters
  names(stopgap.snp.clusters1)[names(stopgap.snp.clusters1)=="start.gwassnp"] <- "clust.init.gwassnp"
  
  stopgap.data.bestld <-   merge(stopgap.data.bestld, stopgap.snp.clusters1,
                                 by = c("snp.ld","snp.gwas","chr.ld", "pos.ld", "pvalue","r2"), all.x = TRUE, all.y=FALSE)
  #sum(is.na(stopgap.data.bestld$cluster))
  
  # Rank the p-values of snp.gwas within each cluster in each pubmedid/disease combination
  stopgap.data.bestld <- stopgap.data.bestld[order(stopgap.data.bestld$pubmedid,
                                                   stopgap.data.bestld$disease,
                                                   stopgap.data.bestld$cluster,
                                                   stopgap.data.bestld$pvalue),]
  stopgap.data.bestld$id <- paste(stopgap.data.bestld$pubmedid, stopgap.data.bestld$disease, stopgap.data.bestld$cluster, sep="_")
  stopgap.data.bestld$id1 <- paste(stopgap.data.bestld$id, stopgap.data.bestld$pvalue, sep="_")
  
  # Add the p-value rank
  tmp <- stopgap.data.bestld[,c("id", "id1")]
  tmp1 <- tmp[match(unique(tmp$id1),tmp$id1),]
  tmp2 <- tmp1[match(unique(tmp1$id),tmp1$id),]

  id.count <- count(tmp1, "id")
  # Put back the order. 
  id.count <- id.count[match(tmp2$id, id.count$id),]
  tmp1$rank <- unlist(lapply(id.count$freq, function(x) 1:x))
  stopgap.data.bestld$gwassnp.p.rank <- tmp1$rank[match(stopgap.data.bestld$id1, tmp1$id1)]
  stopgap.data.bestld$id <- NULL
  stopgap.data.bestld$id1 <- NULL
  rm(tmp)
  rm(tmp1)
  
  save(stopgap.data.bestld, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData", compress = TRUE)
  invisible(stopgap.data.bestld)
  
}


# Generate the bestld data per gene-MeSH.
# Rank by the p-value bins and then by gene scores for  bestld data per gene-MeSH.
# P-value bins (-log10(P)): > 12 as 1; 8-12 as 2 and <8 as 3.
mesh.gene <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData"){
  
  cat("Start subset the bestLD data per each gene and msh combination ...", "\n")
  load(data.file)
  
  stopgap.data.bestld <- subset(stopgap.data.bestld, !is.na(msh))
  # ~ 0.25% rows get removed due to missing msh terms
  
  stopgap.data.bestld$id <- paste(stopgap.data.bestld$gene, stopgap.data.bestld$msh, sep="_")
  
  # p-value bins
  stopgap.data.bestld$p.cat <- cut(stopgap.data.bestld$pvalue, breaks=c(0, 1e-12, 1e-8, 1)) 
  lookup <- data.frame(cat=c("(0,1e-12]", "(1e-12,1e-08]", "(1e-08,1]"), bin=c(1,2,3))
  stopgap.data.bestld$p.bin <- lookup$bin[match(stopgap.data.bestld$p.cat, lookup$cat)]
  
  # Rank the gene_msh by first p-value and then gene score (decreasing)
  stopgap.data.bestld <- stopgap.data.bestld[order(stopgap.data.bestld$id,
                                                   stopgap.data.bestld$p.bin,
                                                   -stopgap.data.bestld$gene.score),]
  #subset(stopgap.data.bestld, id=="A1CF_hyperuricemia")
  
  id.pubmedid.count <- count(stopgap.data.bestld, c("id","pubmedid"))
  id.pubmedid.count <- count(id.pubmedid.count[,-ncol(id.pubmedid.count)], c("id"))
  
  id.cluster.count <- count(stopgap.data.bestld, c("id","cluster"))
  id.cluster.count <- count(id.cluster.count[,-ncol(id.cluster.count)], c("id"))
  
  id.snp.gwas.count <- count(stopgap.data.bestld, c("id","snp.gwas"))
  id.snp.gwas.count <- count(id.snp.gwas.count[,-ncol(id.snp.gwas.count)], c("id"))
  
  id.snp.ld.count <- count(stopgap.data.bestld, c("id","snp.ld"))
  id.snp.ld.count <- count(id.snp.ld.count[,-ncol(id.snp.ld.count)], c("id"))
  
  clus.gene.count <- count(stopgap.data.bestld, c("cluster","gene"))
  clus.gene.count <- count(clus.gene.count[,-ncol(clus.gene.count)], c("cluster"))
  
  stopgap.bestld.gene.mesh <- stopgap.data.bestld[match(unique(stopgap.data.bestld$id), stopgap.data.bestld$id),]
  stopgap.bestld.gene.mesh$asso.count <- id.pubmedid.count$freq[match(stopgap.bestld.gene.mesh$id, id.pubmedid.count$id)]
  stopgap.bestld.gene.mesh$num.cluster <- id.cluster.count$freq[match(stopgap.bestld.gene.mesh$id, id.cluster.count$id)]
  stopgap.bestld.gene.mesh$num.gwas.snp <- id.snp.gwas.count$freq[match(stopgap.bestld.gene.mesh$id, id.snp.gwas.count$id)]
  stopgap.bestld.gene.mesh$num.ld.snp <- id.snp.ld.count$freq[match(stopgap.bestld.gene.mesh$id, id.snp.ld.count$id)]
  stopgap.bestld.gene.mesh$num.gene.cluster <- clus.gene.count$freq[match(stopgap.bestld.gene.mesh$cluster, clus.gene.count$cluster)]
  
  # Add the max gene.score
  # The following is too slow
  #id.max.gene.score <- ddply(stopgap.data.bestld, .(id), summarise, max.gene.score = max(gene.score))
  #stopgap.bestld.gene.mesh$max.gene.score <- id.max.gene.score$max.gene.score[match(stopgap.bestld.gene.mesh$id, id.max.gene.score$id)]  
  # Use the following instead
  id.max.gene.score <- stopgap.data.bestld[,c("id", "gene.score")]
  id.max.gene.score <- id.max.gene.score[order(id.max.gene.score$id, 
                                               -id.max.gene.score$gene.score),]
  id.max.gene.score <- id.max.gene.score[match(unique(id.max.gene.score$id), id.max.gene.score$id),]
  stopgap.bestld.gene.mesh$max.gene.score <- id.max.gene.score$gene.score[match(stopgap.bestld.gene.mesh$id, id.max.gene.score$id)]  
  
  stopgap.bestld.gene.mesh$id <- NULL
  stopgap.bestld.gene.mesh$p.cat <- NULL
  stopgap.bestld.gene.mesh$p.bin <- NULL
  names <- names(stopgap.bestld.gene.mesh)
  names1 <- c("gene", "msh", "msh.tree", "msh.cat", "disease", "pvalue", "gene.score", "gene.rank.max", 
              "gene.rank.min","gwassnp.p.rank", "snp.ld", "snp.gwas", "r2", 
              "cluster", "clust.init.gwassnp")
  stopgap.bestld.gene.mesh <- stopgap.bestld.gene.mesh[,c(names1, names[!(names %in% names1)])]

  system("mkdir Mesh_gene")
  save(stopgap.bestld.gene.mesh, file = "./Mesh_gene/STOPGAP_r2_0.7_bestLD_gene_mesh.RData", compress = TRUE)  
  
  invisible(stopgap.bestld.gene.mesh)
  
}

# Reorder all related files based on columns
# Reorder for the bestLD version data
reorder.data <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData"){
  load(data.file)
  names <- names(stopgap.data.bestld)
  names1 <- c("gene", "msh","msh.tree", "msh.cat", "disease", "pvalue", "gene.score", "gene.rank.max", 
              "gene.rank.min","gwassnp.p.rank", "snp.ld", "snp.gwas", "r2", 
              "cluster", "clust.init.gwassnp")
  stopgap.data.bestld <- stopgap.data.bestld[,c(names1, names[!(names %in% names1)])]
  
  save(stopgap.data.bestld, file = data.file, compress = TRUE) 
  
  invisible(stopgap.data.bestld)
  
}

# Reorder for the full version data
reorder.data1 <- function(data.file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_cluster.RData"){
  load(data.file)
  names <- names(stopgap.data)
  names1 <- c("gene", "msh","msh.tree", "msh.cat", "disease", "pvalue", "gene.score", "gene.rank.max", 
              "gene.rank.min","gwassnp.p.rank", "snp.ld", "snp.gwas", "r2", 
              "cluster", "clust.init.gwassnp")
  stopgap.data <- stopgap.data[,c(names1, names[!(names %in% names1)])]
  
  save(stopgap.data, file = data.file, compress = TRUE) 
  
  invisible(stopgap.data)
  
}

# Rename the R object names by using simpler, descriptive names.
rename.data <- function(){
  
  # gwas
  load("./stopgap_3sources_dbSNP141_msh.RData")
  gwas <- gwas.data
  names(gwas) <- c("snp.gwas.orig", "snp.gwas", "pubmedid","pvalue", "disease", "msh","msh.tree","msh.cat",
                   "init.samp", "source", "rep.samp", "nhlbikey", "loc.paper", "datepub",
                   "gwas.ancestry", "sampsize.dis.rep", "sampsize.dis", "sampsize.rep",
                   "snp.gwasdb","source.gwasdb","hpo.id","hpo.term", "do.id","do.term")
  system("mkdir ../STOPGAP2")
  save(gwas, file = "../STOPGAP2/gwas.RData", compress = TRUE) 
  
  # ld.snps.r2
  load("./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData")
  names(ld.snps.r2) <- c("chr.ld", "pos.ld", "ref.ld", "alt.ld", "chr.gwas", "pos.gwas", "snp.gwas","ref.gwas",
                         "alt.gwas", "snp.ld", "r2", "af.1kg","asn.af.1kg", "amr.af.1kg",
                         "afr.af.1kg", "eur.af.1kg")
  ld.snps.r2 <- ld.snps.r2[,c("snp.ld","chr.ld", "pos.ld", "ref.ld", "alt.ld", "snp.gwas", "chr.gwas", "pos.gwas", "ref.gwas",
                              "alt.gwas", "r2", "af.1kg","asn.af.1kg", "amr.af.1kg",
                              "afr.af.1kg", "eur.af.1kg")]
  save(ld.snps.r2, file = "../STOPGAP2/ld.snps.r2.RData", compress = TRUE) 
  
  # ld.snps
  load("./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData")
  names(ld.snps) <- c("chr.ld", "pos.ld","snp.ld")
  ld.snps <- ld.snps[,c("snp.ld", "chr.ld", "pos.ld")]
  save(ld.snps, file = "../STOPGAP2/ld.snps.RData", compress = TRUE) 
  
  # snp.cluster
  load("./LocusClustering/STOPGAP_SNP_Clusters.RData")
  snp.cluster <- stopgap.snp.clusters[,c("snp.ld","chr.ld","pos.ld","snp.gwas","pvalue","r2","cluster","start.gwassnp")]
  names(snp.cluster)[names(snp.cluster)=="start.gwassnp"] <- "clust.init.gwassnp"
  save(snp.cluster, file = "../STOPGAP2/snp.cluster.RData", compress = TRUE)
  
  # var2gene.simp 
  load("./VarGeneMapping/Var2Gene_vepSimp.RData")
  var2gene.simp <- var2gene.vepSimp
  var2gene.simp <- subset(var2gene.simp, CHR.ld==chr.gencode)
  var2gene.simp$chr.gencode <- NULL
  names(var2gene.simp) <- c("snp.ld", "gene", "chr.ld", "pos.ld", "cat.rdb", "dhs.cor", 
                            "dhs.type","chiapet.cell","chiapet.type","fantom5.tissue","fantom5.type",
                            "eqtl.ref", "eqtl.tissue", "eqtl.scoretype","eqtl.score",
                            "vep.canonical","vep.conseq","vep.aa","vep.condel")
  save(var2gene.simp, file = "../STOPGAP2/var2gene.simp.RData", compress = TRUE)
  
  # stopgap
  load("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_cluster.RData")
  stopgap <- stopgap.data
  #write.table(stopgap,
  #            "../STOPGAP2/stopgap.txt",
  #            sep = "\t", row.names = FALSE, quote = F, na = "")
  save(stopgap, file = "../STOPGAP2/stopgap.RData", compress = TRUE)
  
  # stopgap.bestld
  load("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData")
  stopgap.bestld <- stopgap.data.bestld
  write.table(stopgap.bestld,
              "../STOPGAP2/stopgap.bestld.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")
  save(stopgap.bestld, file = "../STOPGAP2/stopgap.bestld.RData", compress = TRUE) 
  
  # stopgap.gene.mesh
  load("./Mesh_gene/STOPGAP_r2_0.7_bestLD_gene_mesh.RData")
  stopgap.gene.mesh <- stopgap.bestld.gene.mesh
  write.table(stopgap.gene.mesh,
              "../STOPGAP2/stopgap.gene.mesh.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")  
  save(stopgap.gene.mesh, file = "../STOPGAP2/stopgap.gene.mesh.RData", compress = TRUE) 
  
  invisible(stopgap.gene.mesh)
  
}

# Update the snp.cluster data by only keeping three columns "snp.ld", "cluster", "start.gwassnp"
snp.cluster.update <- function(data.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData"){
  load(data.file)
  snp.cluster <- stopgap.snp.clusters[,c("snp.ld", "cluster", "start.gwassnp")]
  names(snp.cluster) <- c("snp.ld", "cluster", "clust.init.gwassnp")
  snp.cluster <- snp.cluster[!duplicated(snp.cluster),]
  save(snp.cluster, file = "../STOPGAP2/snp.cluster.RData", compress = TRUE)
  
  invisible(snp.cluster)
}

# Update the gwas data by adding the chr and pos information
gwas.update <- function(gwas.file = "../STOPGAP2/gwas.RData",
                        rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt"){
  
  load(gwas.file)
  coord <- read.delim(rsid.file,
                      header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  names(coord)[names(coord)=="X.rsID"] <- "rsID"
  
  gwas$GRCh37_Coordinate <- coord$GRCh37_Coordinate[match(gwas$snp.gwas, coord$rsID)]
  gwas$chr.gwas <- unlist(lapply( strsplit(gwas$GRCh37_Coordinate,":"),function(x)x[[1]]))
  gwas$pos.gwas <- unlist(lapply( strsplit(gwas$GRCh37_Coordinate,":"),function(x)x[[2]]))
  col.keep <- c( "snp.gwas.orig", "snp.gwas", "chr.gwas", "pos.gwas", "pubmedid",
                 "pvalue", "disease", "msh","msh.tree","msh.cat","init.samp", "source", "rep.samp", "nhlbikey",
                 "loc.paper", "datepub", "gwas.ancestry", "sampsize.dis.rep", "sampsize.dis",
                 "sampsize.rep", "snp.gwasdb", "source.gwasdb", "hpo.id", "hpo.term", "do.id",
                 "do.term")
  gwas <- gwas[,col.keep]  
  
  save(gwas, file = "../STOPGAP2/gwas.RData", compress = TRUE)
  
  invisible(gwas)
}

# res.update: Update stopgap results 
# Make allele frequencies numeric variables
# Replace "-" with <NA> for VEP annotations
# separate vep.condel into two columns: prediction (character) and score (numeric)
res.update <- function(in.data){

  # Make allele frequencies numeric variables
  for (i in c("af.1kg", "asn.af.1kg", "amr.af.1kg", "afr.af.1kg", "eur.af.1kg")){
    in.data[,i] <- as.numeric(in.data[,i])
  }
  
  # Replace "-" with <NA> for VEP annotations
  # separate vep.condel into two columns: prediction (character) and score (numeric)  
  in.data$vep.canonical[which(in.data$vep.canonical=="-")] <- NA
  in.data$vep.aa[which(in.data$vep.aa=="-")] <- NA
  in.data$vep.condel[which(in.data$vep.condel=="-")] <- NA
  # separate vep.condel into two columns: prediction (character) and score (numeric)
  pred <- unlist(lapply( strsplit(in.data$vep.condel,"\\("),function(x)x[1]))
  score <- unlist(lapply( strsplit(in.data$vep.condel,"\\("),function(x)x[2]))
  score <- as.numeric(gsub("\\)", "", score))
  in.data$vep.condel.pred <- pred
  in.data$vep.condel.score <- score  
  n <- which(names(in.data)=="vep.condel")
  m <- ncol(in.data)
  in.data <- in.data[,c(1:n, (m-1):m , (n+1):(m-2))]  
  in.data$vep.condel <- NULL
  in.data
  
}

# Update all the results for: 
# Make allele frequencies numeric variables
# Replace "-" with <NA> for VEP annotations
# separate vep.condel into two columns: prediction (character) and score (numeric)
results.update <- function(res.path="../STOPGAP2"){
  
  # gwas.RData
  # No need to update

  # var2gene.simp: Replace "-" with <NA> for VEP annotations.
  load(paste(res.path, "var2gene.simp.RData", sep="/"))
  var2gene.simp$vep.canonical[which(var2gene.simp$vep.canonical=="-")] <- NA
  var2gene.simp$vep.aa[which(var2gene.simp$vep.aa=="-")] <- NA
  var2gene.simp$vep.condel[which(var2gene.simp$vep.condel=="-")] <- NA
  # separate vep.condel into two columns: prediction (character) and score (numeric)
  pred <- unlist(lapply( strsplit(var2gene.simp$vep.condel,"\\("),function(x)x[1]))
  score <- unlist(lapply( strsplit(var2gene.simp$vep.condel,"\\("),function(x)x[2]))
  score <- as.numeric(gsub("\\)", "", score))
  var2gene.simp$vep.condel <- NULL
  var2gene.simp$vep.condel.pred <- pred
  var2gene.simp$vep.condel.score <- score
  save(var2gene.simp, file = "../STOPGAP2/var2gene.simp.RData", compress = TRUE)
  
  # stopgap.gene.mesh: 
  # Update (1) mesh terms, (2) Categories
  # Make allele frequencies numeric variables
  # Replace "-" with <NA> for VEP annotations
  # separate vep.condel into two columns: prediction (character) and score (numeric)
  load(paste(res.path, "stopgap.gene.mesh.RData", sep="/"))
  names(stopgap.gene.mesh)[names(stopgap.gene.mesh)=="severity"] <- "vep.severity"
  stopgap.gene.mesh <- res.update(in.data=stopgap.gene.mesh)
  save(stopgap.gene.mesh, file = "../STOPGAP2/stopgap.gene.mesh.RData", compress = TRUE)

  # stopgap.bestld
  load(paste(res.path, "stopgap.bestld.RData", sep="/"))
  names(stopgap.bestld)[names(stopgap.bestld)=="severity"] <- "vep.severity"
  stopgap.bestld <- res.update(in.data=stopgap.bestld)
  write.table(stopgap.bestld,
              "../STOPGAP2/stopgap.bestld.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")  
  save(stopgap.bestld, file = "../STOPGAP2/stopgap.bestld.RData", compress = TRUE)
  
  # stopgap  
  load(paste(res.path, "stopgap.RData", sep="/"))
  names(stopgap)[names(stopgap)=="severity"] <- "vep.severity"
  stopgap <- res.update(in.data=stopgap)
  save(stopgap, file = "../STOPGAP2/stopgap.RData", compress = TRUE)
  
  invisible(stopgap)
  
}


# merge.omim: Merge stopgap.gene.mesh with latest OMIM/Orphanet data processed on 11/13/2104
# Create merge between OMIM, Orphanet, and stopgap.gene.mesh -> stopgap.gwas.omim
# Change "Link" to "pubmedid" and add the OrphID as "pubmedid" when the Source=="Orphanet". 
# Set Rank=1, pval2 as 0 and GeneScore=max(stopgap.gene.mesh$gene.score)
# Remove MSH.TOP and OrphID columms and change the names to being consistent with stopgap.gene.mesh
# Merge the two datasets by using all column names of omim.nomhc data. 
#  Basically it is adding all rows of omim.nomhc to the stopgap.gene.mesh data 
merge.omim <- function(gm.file = "../STOPGAP2/stopgap.gene.mesh.RData",
                       omim.file = "./OMIM_Orphanet/omim.RData"){  
  
  load(gm.file)
  load(omim.file)
  
  omim.1 <- subset(omim, Source=="OMIM")
  omim.2 <- subset(omim, Source=="Orphanet")
  omim.1$Link <- gsub("OMIM:", "", omim.1$Link)
  omim.2$Link <- omim.2$OrphID
  omim <- rbind(omim.1, omim.2)
  omim$MSH.Top <- NULL
  omim$OrphID <- NULL
  omim$Rank <- NULL
  omim$GeneScore <- max(stopgap.gene.mesh$gene.score)
  omim$gene.rank.max <- 1
  omim$gene.rank.min <- 1
  omim$pvalue <- 0
  names(omim) <- c("disease", "pubmedid", "source", "msh", "gene", "gene.score", "gene.rank.max", "gene.rank.min", "pvalue")
  omim <- omim[,c("gene", "msh", "disease", "pvalue", "gene.score", "gene.rank.max", "gene.rank.min", "pubmedid", "source")]
  
  names.orig <- names(stopgap.gene.mesh)
  stopgap.gene.mesh <- merge(stopgap.gene.mesh, omim, 
                             by = names(omim), all.x = TRUE, all.y = TRUE)
  stopgap.gene.mesh <- stopgap.gene.mesh[,names.orig]
  
  write.table(stopgap.gene.mesh,
              "../STOPGAP2/stopgap.gene.mesh.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")  
  save(stopgap.gene.mesh, file = "../STOPGAP2/stopgap.gene.mesh.RData", compress = TRUE)   
  
  invisible(stopgap.gene.mesh)
  
}

# Remove the pubmedid=="22377632" or msh=="transmission" in the stopgap.RData, stopgap.bestld.RData and stopgap.gene.mesh.RData
# #Removing PUBMEDID from GWASDB only: PUBMEDID 21150878 & PUBMEDID 22190364
res.update1 <- function(rm.pubmedid="22377632",
                        rm.pubmedid1=c("21150878","22190364")){
  
  load("../STOPGAP2/stopgap.gene.mesh.RData")
  stopgap.gene.mesh <- subset(stopgap.gene.mesh, pubmedid!=rm.pubmedid)
  stopgap.gene.mesh <- subset(stopgap.gene.mesh, !(pubmedid %in% rm.pubmedid1 & source=="gwasdb") )
  save(stopgap.gene.mesh, file = "../STOPGAP2/stopgap.gene.mesh.RData")
  write.table(stopgap.gene.mesh,
              "../STOPGAP2/stopgap.gene.mesh.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")    
  
  nrow(stopgap.gene.mesh)
  
  load("../STOPGAP2/stopgap.bestld.RData")
  stopgap.bestld <- subset(stopgap.bestld, pubmedid!=rm.pubmedid)
  stopgap.bestld <- subset(stopgap.bestld, !(pubmedid %in% rm.pubmedid1 & source=="gwasdb") )
  save(stopgap.bestld, file = "../STOPGAP2/stopgap.bestld.RData")
  nrow(stopgap.bestld)
  
  load("../STOPGAP2/stopgap.RData")
  stopgap <- subset(stopgap, pubmedid!=rm.pubmedid)
  stopgap <- subset(stopgap, !(pubmedid %in% rm.pubmedid1 & source=="gwasdb") )
  save(stopgap, file = "../STOPGAP2/stopgap.RData")
  nrow(stopgap)
  
}



# Create a function genomics and annotation evidence column for the bestld or gene.mesh data
sg.evid <- function(data) {
  estr <- rep("", nrow(data))
  
  lVep <- !is.na(data$vep.severity)
  lNS <- !is.na(data$vep.condel.pred)
  x <- substr(as.character(data$vep.condel.pred[lNS]), 1, 1)
  estr[lVep] <- paste(estr[lVep], rep("v", sum(lVep)), sep = "")
  estr[lVep & lNS] <- paste(estr[lVep & lNS], x, sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  leQTL <- !is.na(data$eqtl.ref)
  estr[leQTL] <- paste(estr[leQTL], rep("e", sum(leQTL)), sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  lRdb <- !is.na(data$cat.rdb)
  estr[lRdb] <- paste(estr[lRdb], rep("r", sum(lRdb)),
                      substr(as.character(data$cat.rdb)[lRdb], 1, 1),
                      sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  lDhs <- !is.na(data$dhs.cor)
  estr[lDhs] <- paste(estr[lDhs], rep("d", sum(lDhs)), data$dhs.type[lDhs],
                      sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  lChia <- !is.na(data$chiapet.cell)
  estr[lChia] <- paste(estr[lChia], rep("c", sum(lChia)),
                       data$chiapet.type[lChia],
                       sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  lF5 <- !is.na(data$fantom5.tissue)
  estr[lF5] <- paste(estr[lF5], rep("f", sum(lF5)), sep = "")
  
  return(factor(estr))
}

# Add the evidence column to the bestld or gene.mesh data
add.evid <- function(data){
  
  data$evidence <- sg.evid(data)
  
  ## Now merge in the best gene for each association and the evidence
  x <- subset(data, gene.rank.min == 1,
              select = c(gene, disease, snp.gwas, gene.score, evidence))
  x <- subset(x, !duplicated(paste(disease, snp.gwas)))
  names(x) <- paste(names(x), ".best", sep = "")
  data <- merge(data, x,
                by.x = c("disease", "snp.gwas"),
                by.y = c("disease.best", "snp.gwas.best"),
                all = TRUE)
  data
}

# row and col numbers of STOPGAP2 main datasets.
row.col <- function(){
  
  setwd("./STOPGAP2")
  num <- list()
  # gwas
  load("gwas.RData")
  num1 <- c("gwas.RData", nrow(gwas), ncol(gwas))
  num[[1]] <- num1
  
  # ld.snps.r2
  load("ld.snps.r2.RData")
  num2 <- c("ld.snps.r2.RData", nrow(ld.snps.r2), ncol(ld.snps.r2))
  num[[2]] <- num2
  
  # ld.snps
  load("ld.snps.RData")
  num3 <- c("ld.snps.RData", nrow(ld.snps), ncol(ld.snps))
  num[[3]] <- num3
  
  # snp.cluster
  load("snp.cluster.RData")
  num4 <- c("snp.cluster.RData", nrow(snp.cluster), ncol(snp.cluster))
  num[[4]] <- num4
  
  # var2gene.simp 
  load("var2gene.simp.RData")
  num5 <- c("var2gene.simp.RData", nrow(var2gene.simp), ncol(var2gene.simp))
  num[[5]] <- num5
  
  # stopgap
  load("stopgap.RData")
  num6 <- c("stopgap.RData", nrow(stopgap), ncol(stopgap))
  num[[6]] <- num6
  
  # stopgap.bestld
  load("stopgap.bestld.RData")
  num7 <- c("stopgap.bestld.RData", nrow(stopgap.bestld), ncol(stopgap.bestld))
  num[[7]] <- num7
  
  # stopgap.gene.mesh
  load("stopgap.gene.mesh.RData")
  num8 <- c("stopgap.gene.mesh.RData", nrow(stopgap.gene.mesh), ncol(stopgap.gene.mesh))
  num[[8]] <- num8
  
  num <- do.call("rbind", num)
  colnames(num) <- c("file", "Rows", "Cols")
  write.table(num,
              "stopgap2_Rows_Cols.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")
  
  invisible(num)
  
}




