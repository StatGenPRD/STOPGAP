#--------STOPGAP2 R functions ----------------#
# Author: Judong Shen and KIjoung Song

# Trim white spaces
trimWhiteSpace <- function (x) {
  sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

sbind = function(x, y, fill=NA) {
    sbind.fill = function(d, cols){ 
        for(c in cols)
            d[[c]] = fill
        d
    }

    x = sbind.fill(x, setdiff(names(y),names(x)))
    y = sbind.fill(y, setdiff(names(x),names(y)))

    rbind(x, y)
}


# Download GraspFullDataset2.zip data
# url <- c("https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip")
# file <- c("GraspFullDataset2.zip")
# download.file(url, file)

# grasp.import: use unix commands/R scripts to filter the grasp data
# Assume the input data is a *.zip file originally downloaded from the GRASP website.
grasp.import <- function(file = "./Data/GraspFullDataset2.zip", p.thres=1E-04)
{
  cat(paste("Start processing the ", file, " data ...", sep=""), "\n")
  system(paste("zcat", file,  "| head -n 1 > header.txt", sep=" "))
  vars <- names(read.delim("header.txt",
                           header = TRUE, comment.char = "", sep = "\t", as.is = TRUE))
  vars.keep <- c("SNPid.dbSNP134.", "PMID", "Pvalue",
                 "Phenotype","Initial.Sample.Description", "Replication.Sample.Description")
  
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
  if(any(grep("^Gene expression", grasp$Phenotype))) grasp <- grasp[-grep("^Gene expression", grasp$Phenotype),]
  #ignore warning() because it was occurred by '\x'(UTF-8 mode).
  grasp$PMID <- trimWhiteSpace(grasp$PMID)
  save(grasp, file = "grasp.RData", compress = TRUE)
  
  invisible(grasp)
}


# NHGRI.import: filter the NHGRI data
NHGRI.import <- function(file = "./Data/gwas_catalog_v1.0-associations_e85_r2016-08-01.tsv", p.thres=1E-04)
{
  cat(paste("Start processing the NHGRI ", file, " data ...", sep=""), "\n")
  NHGRI <- read.delim(file,
                      header = TRUE, sep = "\t", as.is = TRUE)  
 # NHGRI <- NHGRI[-grep("^[a-z | A-Z]", NHGRI$P.VALUE), ]
  NHGRI <- NHGRI[!is.na(NHGRI$P.VALUE),]
  NHGRI$P.VALUE <- as.numeric(NHGRI$P.VALUE)
  NHGRI <- subset(NHGRI, P.VALUE <= p.thres)
  vars.keep <- c("PUBMEDID","SNPS", "P.VALUE", "DISEASE.TRAIT", "INITIAL.SAMPLE.SIZE", "REPLICATION.SAMPLE.SIZE")
  NHGRI <- NHGRI[,vars.keep]
  names(NHGRI)[names(NHGRI)=="INITIAL.SAMPLE.SIZE"] <- "Initial.Sample.Size"
  names(NHGRI)[names(NHGRI)=="REPLICATION.SAMPLE.SIZE"] <- "Replication.Sample.Size"
  
  # Remove the rows with haplotypes (multiple SNPs) instead of single SNPs
  nhgri <- subset(NHGRI, !grepl(",", SNPS))
  nhgri <- subset(NHGRI, !grepl("x", SNPS))

  save(nhgri, file = "nhgri.RData", compress = TRUE)
  
  invisible(nhgri)
}

# gwasdb.import: filter the gwasdb data
gwasdb.import <- function(file = "./Data/gwasdb_20150819_snp_trait.gz", p.thres=1E-04)
{
  cat(paste("Start processing the gwasdb data ", file, "  ...", sep=""), "\n")
  system(paste("zcat", file,  "| head -n 1 > header.txt", sep=" "))
  vars <- names(read.delim("header.txt",
                           header = TRUE, comment.char = "", sep = "\t", as.is = TRUE))
  vars <- vars[-4]
  gwasdb <- read.delim(gzfile(file),
                       header = FALSE, sep = "\t", as.is = TRUE, skip=1)
  colnames(gwasdb) <- vars 
  #gwasdb <- gwasdb[-grep("^[a-z | A-Z]", gwasdb$P_VALUE), ]
  gwasdb <- gwasdb[!is.na(gwasdb$P_VALUE),]
  gwasdb$P_VALUE <- as.numeric(gwasdb$P_VALUE)  
  
  gwasdb <- subset(gwasdb, P_VALUE <= p.thres)
  if(any(grep("^Gene expression", gwasdb$GWAS_TRAIT))) gwasdb <- gwasdb[-grep("^Gene expression", gwasdb$GWAS_TRAIT),]
  vars.keep <- c("PMID","SNPID.dbSNP","ORI_SNPID", "P_VALUE", "GWAS_TRAIT", "GWAS_INITIAL_SAMPLE_SIZE","SOURCE")
  gwasdb <- gwasdb[,vars.keep]
  gwasdb$SNPID.dbSNP <- trimWhiteSpace(gwasdb$SNPID.dbSNP)
  gwasdb$ORI_SNP <- trimWhiteSpace(gwasdb$ORI_SNPID)
  gwasdb$PMID <- trimWhiteSpace(gwasdb$PMID)
  save(gwasdb, file = "gwasdb.RData", compress = TRUE)
  
  invisible(gwasdb)
}

# phewas.import: use unix commands/R scripts to filter the phewas data
# Assume the input data is a *.csv file originally downloaded from the PheWAS website.
phewas.import <- function(file = "./Data/phewas-catalog.csv", p.thres=1E-04)
{
  cat(paste("Start processing the ", file, " data ...", sep=""), "\n")
  phewas <- read.csv(file,header = TRUE)
  
#  phewas <- phewas[-grep("^[a-z | A-Z]", phewas$p.value), ]
  phewas <- phewas[!is.na(phewas$p.value),]
  phewas$p.value <- as.numeric(phewas$p.value)  
  
  phewas <- subset(phewas, p.value <= p.thres)
  vars.keep <- c("snp","phewas.phenotype","cases", "p.value","phewas.code")
  phewas <- phewas[,vars.keep]
  phewas$snp <- trimWhiteSpace(phewas$snp)
  phewas$phewas.phenotype <- trimWhiteSpace(phewas$phewas.phenotype)
  save(phewas, file = "phewas.RData", compress = TRUE)
  
  invisible(phewas)
}

# venn.plot: venn diagram plot
venn.plot <- function(data.plot,file.name){
#  png(file=file.name,width=10.5,height=10.5,units="in",res=300)
  pdf(file=file.name, height = 8, width = 8)
  par(mar=c(1,1,1,1))
   venn(data.plot)
  dev.off()
}

#binding two data frames that don't have the same set of columns
sbind = function(x, y, fill=NA) {
    sbind.fill = function(d, cols){ 
        for(c in cols)
            d[[c]] = fill
        d
    }

    x = sbind.fill(x, setdiff(names(y),names(x)))
    y = sbind.fill(y, setdiff(names(x),names(y)))

    rbind(x, y)
}

# gwas.filter: filter the gwas data based on disease names and PUBMEDID
gwas.filter <- function(data){
  # Remove the gwas data with some disease names that will not be used 
  no <- grep("methylation levels", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("differential exon", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("differential splicing", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }
  no <- grep("transcript initiation", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("transcript termination", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("gene expression", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }  
  no <- grep("methylation qtl", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("qtl", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("exon skipping", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("differential", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("desialylated glycan peak", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("complex transcript isoform variation", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("recombination", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("expression", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }    
  no <- grep("intron retention", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("transcript levels", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("allele-specific methylation", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("serum metabolite", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }   
  no <- grep("cell lines", trimWhiteSpace(data$disease))
  if (length(no)>0){
    data <- data[-no,]
  }   
 
 x=c("21886157",  # Serum ratio of ...
     "18403759",  # YKL-40 levels
     "21572414",  # "Urinary metabolites"
										 "20037589",  # arg-ptc
										 "22341974", # craniofacial morphology variation 
										 "22595970", # protein abundance levels
										 "19043545",
										 "19798445",
										 "22286219",
                                         "23382691", "23281178","22286219", "17903293", "19043545", "17554300", "19862010","21150878","22190364","22377632")
  # Further remove records in some publications
  data <- subset(data, !(as.character(PUBMEDID) %in% x))
										 
  #   write.table(data,
  #               "stopgap_grasp_nhgri_gwasdb_merged_filtered.txt",
  #               sep = "\t", row.names = FALSE, quote = F, na = "")
  # save(data, file = "stopgap_4sources_clean.RData", compress = TRUE)
  data  
}


# data.merge: merge the grasp, nhgri and gwasdb datasets by first using nhgri records, supplementing 
#              by grasp and gwasdb datasets respectively
# Merge order: NHGRI > GRASP > GWASDB > phewas 
data.merge <- function(nhgri, grasp, gwasdb, phewas){
 cat("Start merging the data from five sources: nhgri, grasp, gwasdb, and PheWAS ...", "\n")
 load("nhgri.RData")
 load("grasp.RData")
 load("gwasdb.RData")
 load("phewas.RData")
 
  nhgri <- rename(nhgri, c(SNPS = 'snp_id', 
                           P.VALUE = "pvalue",
                           DISEASE.TRAIT = 'disease',
                           Initial.Sample.Size = 'Initial.Sample',
                           Replication.Sample.Size = "Replication.Sample"))
  nhgri<-nhgri[,c("snp_id","PUBMEDID","pvalue","disease","Initial.Sample","Replication.Sample")]
  nhgri$pvalue <- as.numeric(format(nhgri$pvalue, digits=3, scientific=T))  
  nhgri <- nhgri[grepl("rs",nhgri$snp_id),]
  nhgri$Source <- "nhgri"
  
  # Remove redundant rows in nhgri data ("PUBMEDID", "snp_id", "pvalue", "disease", "Initial.Sample", "Replication.Sample" "Source" )
  nhgri <- nhgri[!duplicated(nhgri),]
  nhgri <- gwas.filter(nhgri)
  
  grasp <- rename(grasp, c(SNPid.dbSNP134. = 'snp_id', 
                           PMID = 'PUBMEDID',
                           Pvalue = "pvalue",
                           Phenotype = 'disease',
                           Initial.Sample.Description = 'Initial.Sample',
                           Replication.Sample.Description = 'Replication.Sample'))
  grasp<-grasp[,c("snp_id","PUBMEDID","pvalue","disease","Initial.Sample","Replication.Sample")]
  grasp$snp_id <- paste("rs",grasp$snp_id,sep="")
  grasp$Source <- "grasp"
  grasp <- subset(grasp, !is.na(pvalue))
  grasp$pvalue <- as.numeric(format(grasp$pvalue, digits=3, scientific=T)) 
  grasp <- gwas.filter(grasp)
  
  gwasdb <- rename(gwasdb, c(PMID = 'PUBMEDID',
                             SNPID.dbSNP = "gwasdb_SNP_ID",
                             ORI_SNP = 'snp_id', 
                             P_VALUE = "pvalue",
                             GWAS_TRAIT = 'disease',
                             GWAS_INITIAL_SAMPLE_SIZE = 'Initial.Sample',
                             SOURCE = "Source_gwasdb"))
  gwasdb<-gwasdb[,c("snp_id","PUBMEDID","pvalue","disease","Initial.Sample")]
  gwasdb <- gwasdb[grepl("rs",gwasdb$snp_id),]

  # Maually change some of the mis-spelled pvalue here
  gwasdb$pvalue[gwasdb$pvalue=="1.49 E-9"] <- "1.49E-9"
  gwasdb$pvalue[gwasdb$pvalue=="1.05 E-14"] <- "1.05E-14"
  gwasdb$pvalue[gwasdb$pvalue=="1.3E?\\05"] <- "1.3E-05"
  gwasdb$pvalue[gwasdb$pvalue=="1.78E?05"] <- "1.78E-05"
  
  gwasdb <- subset(gwasdb, pvalue!="")
  gwasdb$Source <- "gwasdb" 
  gwasdb$pvalue <- as.numeric(format(gwasdb$pvalue, digits=3, scientific=T)) 
  gwasdb <- gwas.filter(gwasdb)
  
  phewas <- rename(phewas, c(phewas.code= 'PUBMEDID',
                           snp = 'snp_id', 
                           p.value = "pvalue",
                           phewas.phenotype = 'disease',
                           cases='Initial.Sample'))
  phewas <- phewas[,c("snp_id","pvalue","disease","Initial.Sample","PUBMEDID")]
  phewas$PUBMEDID <- paste("PW",phewas$PUBMEDID,sep="")
  phewas <- subset(phewas, !is.na(pvalue))
  phewas <- subset(phewas, !is.na(snp_id))
  phewas <- phewas[grepl("rs",phewas$snp_id),]
  phewas$pvalue <- as.numeric(format(phewas$pvalue, digits=3, scientific=T)) 
  phewas$Source <- "PheWAS"

  grasp$id <- paste(grasp$snp_id, grasp$PUBMEDID,sep="_")
  nhgri$id <- paste(nhgri$snp_id, nhgri$PUBMEDID,sep="_")
  gwasdb$id <- paste(gwasdb$snp_id, gwasdb$PUBMEDID,sep="_")
#  id1 <- unique(grasp$id);  id2 <- unique(nhgri$id); 
#  id3 <- unique(gwasdb$id) 
  
#  data.plot=list(GRASP=id1,NHGRI=id2,GWASDB=id3)
#  file.name="grasp_nhgri_gwasdb_venn.pdf"
#  venn.plot(data.plot,file.name)
  
 # data.plot=list(GRASP=id1,NHGRI=id2,GWASDB=id3,PheWAS=id4,Exome=id5)
 # file.name="grasp_nhgri_gwasdb_phewas_exome_venn.png"
 # venn.plot(data.plot,file.name)
 
  # Subset of grasp data with id (SNP_PMID) not in nhgri data
  grasp1 <- subset(grasp, !(id %in% nhgri$id))
  grasp1$Source <- "grasp"
  # Remove redundant rows in grasp1 data 
  grasp1 <- grasp1[!duplicated(grasp1),]
  
  # Subset of grasp data with id (SNP_PMID) not in nhgri  
  gwasdb1 <- subset(gwasdb, !(id %in% c(nhgri$id, grasp$id)))
  gwasdb1$Source <- "gwasdb"
  # Remove redundant rows in grasp1 data 
  gwasdb1 <- gwasdb1[!duplicated(gwasdb1),]  
 
  data <- sbind(nhgri, grasp1)
  data <- sbind(data, gwasdb1)
  data <- sbind(data, phewas)
  data$id <- paste(data$snp_id, data$PUBMEDID,sep="_")

  data$snp_id <- trimWhiteSpace(data$snp_id)
  data <- data[,-match(c("id"), names(data))]
  gwas.data <- data 
  #grasp has a strange chracter like \U3e32393cs which should convert to "'".
  gwas.data$disease=gsub("\U3e32393cs","'",gwas.data$disease)
  #exclude additional "b" from some snps
  gwas.data$snp_id=gsub("b","",gwas.data$snp_id)
  
  #write.table(data,
  #            "stopgap_grasp_nhgri_gwasdb_merged.txt",
  #            sep = "\t", row.names = FALSE, quote = F, na = "")
#  save(data, file = "stopgap_4sources.RData", compress = TRUE)
  save(gwas.data, file = "stopgap_4sources_clean.RData", compress = TRUE) 

  invisible(data)
  
}

# new.gwassnps: identify the unique GWAS SNPs for LD calculation 
gwas.snps <- function(gwas.file="stopgap_4sources_clean.RData"){ 
  load(gwas.file) 
  snps <- sort(unique(gwas.data$snp_id)) 
  writeLines(snps, "stopgap2_SNPs.txt", sep = "\n", useBytes = FALSE) 
    
  invisible(snps) 
} 

# Coordinate lookup and LD calculation based on the 1KG data for the new SNPs identified in this version
# run.ld: Calculate the LD SNPs based on the 
# rslist is the list of rs IDs with # prefixed comment lines allowed
# https://connect.gsk.com/sites/genetics/GeneticsWIKI/Wiki%20Pages/STOPGAP%20-%20LD%20Calculation.aspx
# before running Get_rsID_coord.py, please add the following path at UNIX
#  PATH=$PATH\:/GWD/bioinfo/projects/GXapp/EntrezDirect/edirect


run.ld <- function(){
  # Input SNPs: stopgap2_SNPs.txt
  
  system("mkdir STOPGAP2_LDResults")
  
  #(1) Coordinate lookup
  system("python2.7 ./Get_rsID_coord.py -o ./STOPGAP_LD -r ./Example/stopgap2_snps.txt") 
    
  #(2) LD Calculation example for chr. 1
  system("./LDcalc.sh 1 ./STOPGAP_LD ./STOPGAP_LD")
  
}

# Update rsIDs in the gwas data (data) to dbSNP141 version by using the ./STOPGAP2_LDResults/rsID_Coordinates.txt file 
# Write out the new GWAS data as "stopgap_4sources_dbSNP141.RData" 
rsID.update <- function(gwas.file="stopgap_4sources.RData",  
                        coor.file="./STOPGAP2_LDResults/rsID_Coordinates.txt") 
                        
{ 
  load(gwas.file) 
  gwas.data <- gwas.filter(data) 
  gwas.data$disease <- tolower(trimWhiteSpace(gwas.data$disease)) 
    
  rsid.coor <- read.delim(coor.file, 
                         header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)  
  names(rsid.coor)[names(rsid.coor)=="X.rsID_searched"] <- "rsID " 
  gwas.data$rsid <- rsid.coor$rsID[match(gwas.data$snp_id, rsid.coor$inputrsIDs )] 
  names <- names(gwas.data)[!(names(gwas.data) %in% c("snp_id", "rsid"))] 
  gwas.data <- gwas.data[,c("snp_id", "rsid", names)] 
    
  save(gwas.data, file = "stopgap_4sources_dbSNP141.RData", compress = TRUE)  
    
  invisible(gwas.data) 
} 


# Find all variants in LD with GWAS variants at r2 > 0.5    
ld.snps <- function(gwas.data="stopgap_4sources_dbSNP141.RData",
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
update.ld.snps.r2 <- function(gwas.file="stopgap_4sources_dbSNP141.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.RData",
                              rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt"){
  cat("Update the ld.snps.r2 data by including the (154526-147706) GWAS SNPs data there ...", "\n")
  load(gwas.file)
  load(ld.file)
  coord <- read.delim(rsid.file,
                      header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  names(coord)[names(coord)=="X.rsIDsearched"] <- "rsID"
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
update1.ld.snps.r2 <- function(gwas.file="stopgap_4sources_dbSNP141.RData",
                              ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.RData",
                              dis = 10000){
  cat("Update the ld.snps.r2 data by removing the LD SNPs more than 10kb away from the GWAS SNPs ...", "\n")
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
                   frq.file="./Data/1KG_AF/gws.frq.RData"){
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

#
#STEP 4. add a genomic evidence resorce to the ld snps
# 
# ld.rdb.dhscor: add the Cat.rdb to the ld snps 
ld.rdb <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                      rdb.file="./Data/genomic_resource/regulomedb14.RData"){
  load(ld.file)
  load(rdb.file)


  # add the Cat.rdb to the ld snps 
  ld.snps$Cat.rdb <- rdb$Cat.rdb[match(ld.snps$SNP.ld, rdb$SNP)]  
  gwas.ld.rdb=ld.snps
  
  save(gwas.ld.rdb, file = "./VarGeneMapping/gwas.ld.rdb.RData", compress = TRUE) 
  
  invisible(gwas.ld.rdb)  
  
}

# For each LD SNP, (1) identify each dhscor where SNP falls within promoter or distal genomic region; 
# (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.dhscor <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                      dhscor.file="./Data/genomic_resource/dhscor.RData",
                      mpos.rng.p = c("Chr.pdhs", "Start.pdhs", "Stop.pdhs"),
                      mpos.rng.d = c("Chr.ddhs", "Start.ddhs", "Stop.ddhs")){
  load(ld.file)
  load(dhscor.file)

  system("mkdir ./VarGeneMapping")  

  dhscor$Chr.ddhs <- gsub("chr", "", dhscor$Chr.ddhs)
  dhscor$Chr.pdhs <- gsub("chr", "", dhscor$Chr.pdhs)
  dhscor$Gene <- as.character(dhscor$Gene)
  dhscor.p <- dhscor[,c("Chr.pdhs", "Start.pdhs", "Stop.pdhs", "Gene")]
  dhscor.d <- dhscor[,c("Chr.ddhs", "Start.ddhs", "Stop.ddhs", "Gene", "fpr")]
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    # promoter
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    dhscor.p.i <- subset(dhscor.p, Chr.pdhs==i)
    dhscor.p.i$fpr <- 0
    
    ## Expand range as one row per position for dhscor.p.i
    dhscor.p.i$Length <- dhscor.p.i[, mpos.rng.p[3]] - dhscor.p.i[, mpos.rng.p[2]] + 1
    dhscor.p.i.exp <- data.frame(Position = as.vector(unlist(apply(dhscor.p.i[, mpos.rng.p[2:3]], 1,
                                                                   function(x) x[1]:x[2]))),
                                 Gene = rep(as.character(dhscor.p.i$Gene), dhscor.p.i$Length),
                                 fpr = rep(dhscor.p.i$fpr, dhscor.p.i$Length))
    dhscor.p.i.exp <- dhscor.p.i.exp[!duplicated(dhscor.p.i.exp),]
    tmp <- merge(ld.snps.i, dhscor.p.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp1 <- tmp[!duplicated(tmp),]
    tmp1$type <- NA
    tmp1$type[!is.na(tmp1$fpr)] <- "p"
    
    # distal
    dhscor.d.i <- subset(dhscor.d, Chr.ddhs==i)
    
    ## Expand range as one row per position for dhscor.d.i
    dhscor.d.i$Length <- dhscor.d.i[, mpos.rng.d[3]] - dhscor.d.i[, mpos.rng.d[2]] + 1
    dhscor.d.i.exp <- data.frame(Position = as.vector(unlist(apply(dhscor.d.i[, mpos.rng.d[2:3]], 1,
                                                                   function(x) x[1]:x[2]))),
                                 Gene = rep(as.character(dhscor.d.i$Gene), dhscor.d.i$Length),
                                 fpr = rep(dhscor.d.i$fpr, dhscor.d.i$Length))
    dhscor.d.i.exp <- dhscor.d.i.exp[!duplicated(dhscor.d.i.exp),]
    tmp <- merge(ld.snps.i, dhscor.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$fpr)] <- "d"
    
    tmp <- merge(tmp1, tmp2, 
                 by = c("CHR.ld","POS.ld", "SNP.ld","Gene",  "fpr", "type"), 
                 all.x = TRUE, all.y = TRUE)
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  res$id=paste(res$SNP.ld,res$Gene,res$type,sep="_")

  c4<-aggregate(fpr ~ id, res,mean)
  res=res[!duplicated(res$id),]
  res$fpr=NULL
  res=merge(res,c4,by="id",all.x=T,all.y=F)
  res$id <- NULL
  gwas.ld.dhscor <- res
#  system("mkdir VarGeneMapping")
  save(gwas.ld.dhscor, file = "./VarGeneMapping/gwas.ld.dhscor.RData", compress = TRUE) 
  
  invisible(res)  
  
}

# ----- Merge the CHIA-PET data
# ld.chiapet: For each LD SNP, (1) indentify each chiapet where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.chiapet <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                       chiapet.file="./Data/genomic_resource/CHIA.PET_5celltypes.RData",
                       mpos.rng.d = c("dChr", "dStart", "dEnd")){
  cat("Start merging CHIA-PET data...", "\n")
  load(ld.file)
  load(chiapet.file)
  names(CHIA.PET)[names(CHIA.PET)=="genes"] <- "Gene"
  CHIA.PET$CellType <- gsub(".bed", "", CHIA.PET$CellType)
  CHIA.PET$Gene <- as.character(CHIA.PET$Gene)
  CHIA.PET.d <- CHIA.PET[,c("dChr", "dStart", "dEnd", "Gene", "CellType")]
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
 
    # distal
    CHIA.PET.d.i <- subset(CHIA.PET.d, dChr==i)
    ld.snps.i <- subset(ld.snps, CHR.ld==i)

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
    
    tmp <- tmp2
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  
  gwas.ld.chiapet <- res
  save(gwas.ld.chiapet, file = "./VarGeneMapping/gwas.ld.chiapet.RData", compress = TRUE) 
  
  invisible(res)  
  
}


# ld.eQTL: For each LD SNP, (1) identify each eQTL where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.eQTL <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                    eQTL.file="./Data/genomic_resource/eQTL_single_44Tissues.RData"){
  cat("Start merging eQTL data...", "\n")
  load(ld.file)
  load(eQTL.file)

 # "rs_id_dbSNP142_GRCh37p13","gene_name","gene_type","Tissue","Score"
  eQTL <- eQTL[,c("ID", "hg19", "tissue.pvalue")]
  names(eQTL) <- c("SNP", "Gene", "Tissue.abb.pvalue")
 
  gwas.ld.eQTL.s <- merge(ld.snps, eQTL,  by.x = "SNP.ld",
                        by.y = "SNP", all.x = TRUE, all.y = FALSE)
  
  save(gwas.ld.eQTL.s, file = "./VarGeneMapping/gwas.ld.eQTL.single.RData", compress = TRUE) 
  
  invisible(gwas.ld.eQTL.s)  
  
}

ld.multi.eQTL <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData" ,
                          eQTL.multi.file="./Data/genomic_resource/eQTL_multi_9Tissues.RData"){

  cat("Start merging eQTL data...", "\n")
  load(ld.file)
#   Gene chrom           Gene.ID    start      end strand length  source
  load(eQTL.multi.file)
     
  gwas.ld.eQTL.m <- merge(ld.snps, meQTL,  by.x = "SNP.ld",
                        by.y = "snp", all.x = TRUE, all.y = FALSE)
#  gwas.ld.eQTL.m$Tissue <- ifelse(is.na(gwas.ld.eQTL.m$Gene),NA,"AS;AT;HLV;Lu;MS;NT;SSEL;Th;WB")
  gwas.ld.eQTL.m=gwas.ld.eQTL.m[,c("SNP.ld","CHR.ld","POS.ld","Gene")]
  save(gwas.ld.eQTL.m, file = "./VarGeneMapping/gwas.ld.eQTL.multi.RData", compress = TRUE) 
  
  invisible(gwas.ld.eQTL.m)  
  
}


# ld.fantom5: For each LD SNP, (1) indentify each fantom5 where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.fantom5 <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",  
                       fantom5.file="./Data/genomic_resource/FANTOM5.RData",
                       mpos.rng.d = c("chr", "startpos", "endpos")){
  load(ld.file)
  load(fantom5.file)
  names(fantom5)[names(fantom5)=="gene"] <- "Gene"
  names(fantom5)[names(fantom5)=="X.chrom"] <- "chr"
  fantom5$Gene <- as.character(fantom5$Gene)
  fantom5.d <- fantom5[,c("chr", "startpos", "endpos", "Gene", "Transcript", "fpr")]
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
                                  Tissue = rep(fantom5.d.i$Transcript, fantom5.d.i$Length),
								  fpr = rep(fantom5.d.i$fpr, fantom5.d.i$Length))
    fantom5.d.i.exp <- fantom5.d.i.exp[!duplicated(fantom5.d.i.exp),]
    tmp <- merge(ld.snps.i, fantom5.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$Tissue)] <- "d"
    tmp <- tmp2[,c("CHR.ld","POS.ld", "SNP.ld", "Gene",  "Tissue", "fpr", "type")] 
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  gwas.ld.fantom5 <- res
  save(gwas.ld.fantom5, file = "./VarGeneMapping/gwas.ld.fantom5.RData", compress = TRUE) 
  
  invisible(res)  
    
}
# ld.CHIC: For each LD SNP, (1) indentify each CHIC where SNP falls within promoter or distal genomic region; 
#                           (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.CHIC <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                       CHIC.file="./Data/genomic_resource/chic.RData",
                       mpos.rng.d = c("chr", "startpos", "endpos")){
  load(ld.file)
  load(CHIC.file)
 
  CHIC$Gene <- as.character(CHIC$Gene)
  CHIC.d <- CHIC[,c("chr", "startpos", "endpos", "Gene", "tissue")]
  CHIC.d$chr <- gsub("chr", "", CHIC.d$chr)
  
  ## Run it by chromosome to make it easier memory-wise
  chr <- unique(ld.snps$CHR.ld)
  
  cpos<- list()
  for(i in chr) {
    cat(i, "\n")
    # no promoter
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    
    # distal
    CHIC.d.i <- subset(CHIC.d, chr==i)
    
    ## Expand range as one row per position for CHIC.d.i
    CHIC.d.i$Length <- CHIC.d.i[, mpos.rng.d[3]] - CHIC.d.i[, mpos.rng.d[2]] + 1
    CHIC.d.i.exp <- data.frame(Position = as.vector(unlist(apply(CHIC.d.i[, mpos.rng.d[2:3]], 1,
                                                                    function(x) x[1]:x[2]))),
                                  Gene = rep(as.character(CHIC.d.i$Gene), CHIC.d.i$Length),
                                  Tissue = rep(CHIC.d.i$tissue, CHIC.d.i$Length))
    CHIC.d.i.exp <- CHIC.d.i.exp[!duplicated(CHIC.d.i.exp),]
    tmp <- merge(ld.snps.i, CHIC.d.i.exp, by.x = "POS.ld",
                 by.y = "Position", all.x = TRUE, all.y = FALSE)
    tmp2 <- tmp[!duplicated(tmp),]
    tmp2$type <- NA
    tmp2$type[!is.na(tmp2$Tissue)] <- "d"
    tmp <- tmp2[,c("CHR.ld","POS.ld", "SNP.ld", "Gene",  "Tissue", "type")] 
    
    cpos[[i]] <- tmp
  }
  res <- do.call("rbind", cpos)
  gwas.ld.CHIC <- res
  save(gwas.ld.CHIC, file = "./VarGeneMapping/gwas.ld.CHIC.RData", compress = TRUE) 
  
  invisible(res)  
  
}

# add cato to ld snps 
ld.cato <-function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                   cato.file="./Data/genomic_resource/cato.RData"){
				   
    load(ld.file)
	load(cato.file)
	
# "chrom"           "snpChromStart"   "snpChromEnd"     "rsid"
# "pred.fit.pctSig" "strand"          "motifname"       "position"
# "ref.allele"      "nonref.allele"   "Cell_types"

  cato <- cato[,c("chrom","rsid","pred.fit.pctSig","motifname","Cell_types")]

  gwas.ld.cato <- merge(ld.snps, cato,  by.x = "SNP.ld",
                        by.y = "rsid", all.x = TRUE, all.y = FALSE)
  
  save(gwas.ld.cato, file = "./VarGeneMapping/gwas.ld.cato.RData", compress = TRUE) 
  
  invisible(gwas.ld.cato)  

}

# gwas.ld.deltaSVM.RData is generated and downloaded at Data/genomic_resource 
ld.deltaSVM <-function(){

  system("cp ./Data/genomic_resource/gwas.ld.deltaSVM.RData ./VarGeneMapping")

}

# gwas.ld.phylop.RData is generated and downloaded at Data/genomic_resource 
ld.phylop <-function(){

  system("cp ./Data/genomic_resource/gwas.ld.phylop.RData ./VarGeneMapping")

}


## nearest.gene: Take a merged data set of genes and SNPs on a single
##  chromosome, sorted by position (optional; this gets repeated) with
##  an indicator for whether the entry is a SNP or a Gene, and identify
##  the nearest gene.
nearest.gene <- function(snp.gene)
{
  ## Make sure it is sorted by position
  snp.gene <- snp.gene[order(snp.gene$POS),]
  snps <- grep("SNP", snp.gene$type)
  genes <- grep("Gene", snp.gene$type)
  genes <- c(genes, snps[length(snps)] + 1) # Add dummy gene to end

  ## Initialize previous and next gene values for first variant
  prev.gene <- next.gene <- numeric(length(snps))
  if(genes[1] > snps[1]){
      prev.gene[1] <- NA
      x <- 0
  }else {
      prev.gene[1] <- genes[x <- which(genes == max(genes[genes < snps[1]],
                                           na.rm = TRUE))]
  }
  next.gene[1] <- genes[x + 1]

  for(i in 2:length(snps)) {
    if(next.gene[i - 1] > snps[i]) { # Same genes as for previous SNP
      prev.gene[i] <- prev.gene[i - 1]
      next.gene[i] <- next.gene[i - 1]
    } else {  # Figure out where the new SNP position falls
     prev.gene[i] <- genes[x <- which(genes == max(genes[genes < snps[i]],
                                           na.rm = TRUE))]
      next.gene[i] <- genes[x + 1]
    }
  }

  ## Remove the trailing dummy gene
  genes <- genes[-length(genes)]

  ## Calculate distances to adjacent genes
  prev.dist <- snp.gene$POS[snps] - snp.gene$POS[prev.gene]
  next.dist <- snp.gene$POS[next.gene] - snp.gene$POS[snps]

  snp.gene$prev.gene <- snp.gene$prev.dist <- rep(NA, nrow(snp.gene))
  snp.gene$prev.gene[snps] <- snp.gene$ID[prev.gene]
  snp.gene$prev.dist[snps] <- prev.dist

  snp.gene$next.gene <- snp.gene$next.dist <- rep(NA, nrow(snp.gene))
  snp.gene$next.gene[snps] <- snp.gene$ID[next.gene]
  snp.gene$next.dist[snps] <- next.dist

  ## Keep only SNPs
  snp.gene <- snp.gene[snps,]

  ## Get nearest gene and distance
  near.ind <- apply(snp.gene[, c("prev.dist", "next.dist")], 1,
                    function(x) which(x == min(x, na.rm = TRUE)))
  near.ind <- as.integer(gsub("1:2","1",near.ind))
  snp.gene$near.dist <- snp.gene$near.gene <- rep(NA, nrow(snp.gene))
  snp.gene$near.gene[near.ind == 1] <- snp.gene$prev.gene[near.ind == 1]
  snp.gene$near.gene[near.ind == 2] <- snp.gene$next.gene[near.ind == 2]
  snp.gene$near.dist[near.ind == 1] <- snp.gene$prev.dist[near.ind == 1]
  snp.gene$near.dist[near.ind == 2] <- snp.gene$next.dist[near.ind == 2]

  invisible(snp.gene)
}

#fing gene location based on gencode V19
ld.posgene <- function(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                       Posgene.file="/Data/genomic_resource/gencode_v19_clean.RData"){
  load(ld.file)
  load(Posgene.file) # source: gencode=20224 
  gencode$Gene=trimWhiteSpace(gencode$Gene)

  ch.d=c(as.character(1:22),"X")

  all.result=list()
  k=0
for(i in ch.d){ 
  k<-k+1
    ld.snps.i <- subset(ld.snps, CHR.ld==i)
    ld.snps.i <- ld.snps.i[order(ld.snps.i$POS.ld),]

    gencode.i <- subset(gencode, chrom==i)
    gencode.i <- gencode.i[order(gencode.i$start),]

	x1=ld.snps.i[,c("POS.ld","SNP.ld")]
	x1$type="SNP"
	names(x1)=c("POS","ID","type")
	x2=gencode.i[,c("Gene","start")]
	x2$type="Gene"
	names(x2)=c("ID","POS","type")
	tmp=rbind(x1,x2)
 
    tmp=nearest.gene(tmp)
 
   tmp <- tmp[,c("ID","POS","near.gene")]
   names(tmp) <- c("SNP.ld","POS.ld","Gene")
   tmp$CHR.ld=i   
  all.result[[k]] <- tmp
  }
  gwas.ld.posgene <- do.call("rbind", all.result)
  save(gwas.ld.posgene, file ="./VarGeneMapping/gwas.ld.posgene.RData", compress = TRUE) 

  invisible(gwas.ld.posgene)  
}

# ld.vep: Merge gwas ld SNPs and VEP information. 
# gwas.ld.vep.RData is generated and downloaded at Data/genomic_resource 
ld.vep <-function(){

  system("cp ./Data/genomic_resource/gwas.ld.vep.full.RData ./VarGeneMapping")

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
  gwas.ld.vep3$severity <- severity$Severity[match(gwas.ld.vep3$Consequence, severity$Consequences)]
  gwas.ld.vep3 <- gwas.ld.vep3[order(gwas.ld.vep3$id, gwas.ld.vep3$severity),]
  gwas.ld.vep3.1 <- gwas.ld.vep3[match(unique(gwas.ld.vep3$id), gwas.ld.vep3$id),]
  
  gwas.ld.vep4 <- subset(gwas.ld.vep2, !(CANONICAL %in% "YES"))
  gwas.ld.vep5 <- subset(gwas.ld.vep4, !(id %in%  gwas.ld.vep3.1$id))
  gwas.ld.vep5$severity <- severity$Severity[match(gwas.ld.vep5$Consequence, severity$Consequences)]
  gwas.ld.vep5 <- gwas.ld.vep5[order(gwas.ld.vep5$id, gwas.ld.vep5$severity),]
  gwas.ld.vep6 <- gwas.ld.vep5[match(unique(gwas.ld.vep5$id), gwas.ld.vep5$id),]
  
  gwas.ld.vep1$id <- NULL
  gwas.ld.vep3.1$id <- NULL
  gwas.ld.vep3.1$severity <- NULL
  gwas.ld.vep6$severity <- NULL
  gwas.ld.vep6$id <- NULL
  
  gwas.ld.vep.simplify <- rbind(gwas.ld.vep1, gwas.ld.vep3.1, gwas.ld.vep6)
  # sum(duplicated(gwas.ld.vep.simplify))  # 0
  save(gwas.ld.vep.simplify, file = "./VarGeneMapping/gwas.ld.vep.simplify.RData", compress = TRUE) 
  
  invisible(gwas.ld.vep.simplify)  
}

#merge.ld.simp: Merge based on SNP.ld and gene from each data
# % of unique (ld) GWAS SNPs with at least one genes
merge.ld.simp <- function(dhscor.file ="./VarGeneMapping/gwas.ld.dhscor.RData",
                     chiapet.file = "./VarGeneMapping/gwas.ld.chiapet.RData",
                     fantom5.file = "./VarGeneMapping/gwas.ld.fantom5.RData",
					 CHIC.file="./VarGeneMapping/gwas.ld.CHIC.RData",
                     eQTL.s.file = "./VarGeneMapping/gwas.ld.eQTL.single.RData",
                     eQTL.m.file = "./VarGeneMapping/gwas.ld.eQTL.multi.RData",
					 cato.file="./VarGeneMapping/gwas.ld.cato.RData",
					 delta.file="./VarGeneMapping/gwas.ld.deltaSVM.RData",
					 phylop.file="./VarGeneMapping/gwas.ld.phylop.RData",
                     posgene.file = "./VarGeneMapping/gwas.ld.posgene.RData",
                     vep.file.simplified = "./VarGeneMapping/gwas.ld.vep.simplify.RData",
				     Posigene.file="./Data/gencode_v19_clean.RData")
  cat("Start merging the dhscor, chip-pet, fantom5 and vep data ...", "\n")
  load(dhscor.file)
  load(chiapet.file)
  load(fantom5.file)
  load(CHIC.file)
  load(eQTL.s.file)
  load(eQTL.m.file)
  load(cato.file)
  load(delta.file)
  load(phylop.file)
  load(posgene.file)
  load(vep.file.simplified)
  load(Posigene.file)
  
#  gwas.ld.rdb <- subset(gwas.ld.rdb,!duplicated(SNP.ld))
  gwas.ld.dhscor <- rename(gwas.ld.dhscor, c(fpr = 'DHS.fpr',
                                             type = "DHS.Type"))
   gwas.ld.chiapet <- rename(gwas.ld.chiapet, c(CellType = 'CHIAPET.CellType',
                                               type = "CHIAPET.Type"))
  gwas.ld.fantom5 <- rename(gwas.ld.fantom5, c(Tissue = 'FANTOM5.Tissue',
                                               fpr = 'FANTOM5.fpr',
                                               type = "FANTOM5.Type"))
  gwas.ld.CHIC <- rename(gwas.ld.CHIC, c(Tissue = 'CHIC.Tissue',
                                           type = "CHIC.Type"))
  gwas.ld.eQTL.s <- rename(gwas.ld.eQTL.s, c(Tissue.abb.pvalue = "Single.eQTL.Tissue.p_value"))
  gwas.ld.eQTL.m <- rename(gwas.ld.eQTL.m, c(Tissue = 'Multi.eQTL.Tissue'))
  gwas.ld.cato <- rename(gwas.ld.cato, c(pred.fit.pctSig ="cato.score",
                                         motifname ="cato.motifname",
										 Cell_types="cato.CellType"))
  gwas.ld.deltaSVM <- rename(gwas.ld.deltaSVM, c(top5="deltaSVM.Top5.CellType",
                                                 q95="deltaSVM.quantile.95p"))
  gwas.ld.phylop=gwas.ld.phylop[,c("SNP.ld","POS.ld","CHR.ld","score")]														 
  gwas.ld.phylop <- rename(gwas.ld.phylop, c(score="phyloP.score"))
  gwas.ld.posgene <- rename(gwas.ld.posgene, c(Gene = 'Gene')) 
  gwas.ld.vep.simplify <- rename(gwas.ld.vep.simplify, c(SYMBOL = 'Gene', 
                                                         CANONICAL = 'VEP.Canonical',
                                                         Consequence = "VEP.Consequence",
                                                         Amino_acids = "VEP.AA",
                                                         Condel = "VEP.Condel"))

  gwas.ld.chiapet$CHIAPET.CellType <-trimWhiteSpace(gwas.ld.chiapet$CHIAPET.CellType)  
  gwas.ld.chiapet$CHIAPET.Type <-trimWhiteSpace(gwas.ld.chiapet$CHIAPET.Type)
  gwas.ld.eQTL.s$Single.eQTL.Tissue.p_value <- trimWhiteSpace(gwas.ld.eQTL.s$Single.eQTL.Tissue.p_value)
  gwas.ld.eQTL.s$Single.eQTL.GeneType <- trimWhiteSpace(gwas.ld.eQTL.s$Single.eQTL.GeneType)
  
  gencode=gencode[,c("Gene","start","end")]

  gwas.ld.dhscor1 <- subset(gwas.ld.dhscor, (!is.na(Gene) & !is.na(DHS.fpr)))
  gwas.ld.dhscor1$Gene <- trimWhiteSpace(gwas.ld.dhscor1$Gene)
  gwas.ld.dhscor1$SNP.ld <- trimWhiteSpace(gwas.ld.dhscor1$SNP.ld)

  gwas.ld.fantom51 <- subset(gwas.ld.fantom5, (!is.na(Gene) & !is.na(FANTOM5.fpr)))
  gwas.ld.fantom51$Gene <- trimWhiteSpace(gwas.ld.fantom51$Gene)
  gwas.ld.fantom51$SNP.ld <- trimWhiteSpace(gwas.ld.fantom51$SNP.ld)

  gwas.ld.chiapet1 <- subset(gwas.ld.chiapet, !is.na(Gene))
  gwas.ld.chiapet1$id<-paste(gwas.ld.chiapet1$Gene,gwas.ld.chiapet1$SNP.ld,sep="_")
  id.count <- count(gwas.ld.chiapet1, "id")
  colnames(id.count)<-c("id","CHIAPET.freq")
  c1<-aggregate(CHIAPET.CellType ~ id, gwas.ld.chiapet1,function(x){y=table(x); names(y)[y==max(y)]}) 
  c1$CHIAPET.CellType=gsub("\"","",c1$CHIAPET.CellType)
  c1$CHIAPET.CellType=gsub("c","",c1$CHIAPET.CellType)
  c1$CHIAPET.CellType=gsub(", ",";",c1$CHIAPET.CellType)
  c12<-merge(c1,id.count,by="id")
  gwas.ld.chiapet1<-gwas.ld.chiapet1[!duplicated(gwas.ld.chiapet1$id),]
  gwas.ld.chiapet1$CHIAPET.CellType<-NULL
  gwas.ld.chiapet1<-merge(gwas.ld.chiapet1,c12,by="id")
  gwas.ld.chiapet1$id<-NULL
  gwas.ld.chiapet1$Gene <- trimWhiteSpace(gwas.ld.chiapet1$Gene)
  gwas.ld.chiapet1$SNP.ld <- trimWhiteSpace(gwas.ld.chiapet1$SNP.ld)

  gwas.ld.CHIC1 <- subset(gwas.ld.CHIC, !is.na(Gene))
  gwas.ld.CHIC1$Gene <- trimWhiteSpace(gwas.ld.CHIC1$Gene)
  gwas.ld.CHIC1$SNP.ld <- trimWhiteSpace(gwas.ld.CHIC1$SNP.ld)

  gwas.ld.eQTL.s1 <- subset(gwas.ld.eQTL.s, !is.na(Gene))
  gwas.ld.eQTL.s1$Gene <- trimWhiteSpace(gwas.ld.eQTL.s1$Gene)
  gwas.ld.eQTL.s1$SNP.ld <- trimWhiteSpace(gwas.ld.eQTL.s1$SNP.ld)

  gwas.ld.eQTL.m1 <- subset(gwas.ld.eQTL.m, !is.na(Gene))
  gwas.ld.eQTL.m1$Gene <- trimWhiteSpace(gwas.ld.eQTL.m1$Gene)
  gwas.ld.eQTL.m1$SNP.ld <- trimWhiteSpace(gwas.ld.eQTL.m1$SNP.ld)
  gwas.ld.eQTL.m1$Tissue <-"multi_eQTL"

  gwas.ld.posgene1 <- subset(gwas.ld.posgene, !is.na(Gene))
  gwas.ld.posgene1$Gene <- trimWhiteSpace(gwas.ld.posgene1$Gene)
  gwas.ld.posgene1$SNP.ld <- trimWhiteSpace(gwas.ld.posgene1$SNP.ld)

  gwas.ld.vep.simplify1 <- subset(gwas.ld.vep.simplify, !is.na(Gene) & gwas.ld.vep.simplify$Gene!="-")
  gwas.ld.vep.simplify1$Gene <- trimWhiteSpace(gwas.ld.vep.simplify1$Gene)
  gwas.ld.vep.simplify1$SNP.ld <- trimWhiteSpace(gwas.ld.vep.simplify1$SNP.ld)

  gwas.ld.rdb$SNP.ld <- trimWhiteSpace(gwas.ld.rdb$SNP.ld)
  gwas.ld.cato$SNP.ld <- trimWhiteSpace(gwas.ld.cato$SNP.ld)
  gwas.ld.deltaSVM$SNP.ld <- trimWhiteSpace(gwas.ld.deltaSVM$SNP.ld)
  gwas.ld.phylop$SNP.ld <- trimWhiteSpace(gwas.ld.phylop$SNP.ld)
  gencode$Gene <- trimWhiteSpace(gencode$Gene)

#  gwas.ld.vep.simplify1 <- subset(gwas.ld.vep.simplify1, (!is.na(VEP.Condel) & gwas.ld.vep.simplify1$VEP.Condel !="_") | (!is.na(VEP.Consequence) & gwas.ld.vep.simplify1$VEP.Consequence !="_"))
# rdb, cato, deltaSVM and phylop have No Gene information!!!
 
  gwas.ld.d <- merge(gwas.ld.dhscor1, gwas.ld.chiapet1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.fantom51, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.CHIC1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.eQTL.s1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.eQTL.m1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.vep.simplify1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.posgene1, 
                     by = c("SNP.ld", "Gene","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.rdb, 
                     by = c("SNP.ld","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.cato, 
                     by = c("SNP.ld","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.deltaSVM, 
                     by = c("SNP.ld","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gwas.ld.phylop, 
                     by = c("SNP.ld","CHR.ld", "POS.ld"), all = TRUE)
  gwas.ld.d <- merge(gwas.ld.d, gencode, by="Gene",all.x=TRUE,all.y=FALSE)
  
  gwas.ld.d$distance <- ifelse(gwas.ld.d$POS.ld >= gwas.ld.d$start & gwas.ld.d$POS.ld <= gwas.ld.d$end,0,ifelse(gwas.ld.d$POS.ld >gwas.ld.d$start & gwas.ld.d$POS.ld > gwas.ld.d$end, gwas.ld.d$POS.ld-gwas.ld.d$end,gwas.ld.d$start-gwas.ld.d$POS.ld))


# Add the CHR.ld, POS.ld back
#  gwas.ld.d$CHR.ld <- gwas.ld.rdb.dhscor$CHR.ld[match(gwas.ld.d$SNP.ld, gwas.ld.rdb.dhscor$SNP.ld)]
#  gwas.ld.d$POS.ld <- gwas.ld.rdb.dhscor$POS.ld[match(gwas.ld.d$SNP.ld, gwas.ld.rdb.dhscor$SNP.ld)]
  
    gwas.ld.d <- rename(gwas.ld.d, c(start = 'Gene.start',
	                                 end ="Gene.end", 
                                     distance = "distance.from.gene"))

   var2gene.simplified <- gwas.ld.d

  save(var2gene.simplified, file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData", compress = TRUE)
 

  invisible(var2gene.simplified)
}

#We restricted our mappings to protein coding genes listed in GENCODE or RefSeq.
# Subset of the var2gene data by using the union of RefSeq and Gencode gene names
# Also Update the var2gene mapping data by removing any rows with conflicting chr.ld with those from the chr information from gencode data
# 
var2gene.filter <- function(var2gene.file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData",
                            gencode.file = "./Data/gencode_v19_clean.RData"){
  
  cat("Start to filter var2gene data ...", "\n")
  load(var2gene.file)
  load(gencode.file)
  
  # Filter down the gene names not in Gencode gene names
  # exclude all non-protein coding genes from var2gene.simplified data
  var2gene.vepSimp <- subset(var2gene.simplified, (trimWhiteSpace(Gene) %in% trimWhiteSpace(gencode$Gene)) & (CHR.ld %in% gencode$chrom))
  var2gene.vepSimp <- subset(var2gene.vepSimp, !is.na(Gene))  
  
  save(var2gene.vepSimp, file = "./VarGeneMapping/Var2Gene_vepSimp.RData", compress = TRUE)
  #save(var2gene.vepFull, file = "./VarGeneMapping/Var2Gene_vepFull.RData", compress = TRUE)
  invisible(var2gene.vepSimp)
}

# gwas.mesh: add mesh terms to the GWAS data. 
# Run it on Windows only
gwas.mesh <- function(gwas.file="./stopgap_4sources_dbSNP141.RData",
                      disease2mesh.file = "./Data/MeSH_Aug2016.txt"){
  
  cat("Start adding mesh terms to the GWAS data ...", "\n")
  load(gwas.file)
  dis2gene <- read.delim(disease2mesh.file,
                         header = TRUE, comment.char = "", na.strings = c("",".",NA), sep = "\t", as.is = TRUE) 
  
  gwas.data$disease <- iconv(gwas.data$disease,"WINDOWS-1252","UTF-8")
  gwas.data$disease <- tolower(trimWhiteSpace(gwas.data$disease))  
#  gwas.data <- gwas.filter(gwas.data)
  dis2gene$disease <- iconv(dis2gene$disease,"WINDOWS-1252","UTF-8")
  dis2gene$disease <- tolower(dis2gene$disease)
  dis2gene$disease <- trimWhiteSpace(dis2gene$disease)
  
  #  all(unique(gwas.data$disease) %in% dis2gene$disease)  # FALSE, two disease names look so weird, remove them. 
  # unique(gwas.data$disease)[!(unique(gwas.data$disease) %in% unique(dis2gene$disease))]
  gwas.data <- subset(gwas.data, disease %in%  unique(dis2gene$disease))
  # on Unix, the tolower() will have error: Error in tolower(gwas.data$disease) : invalid multibyte string 3997
  # Run on Windows
  
  dis2gene$MeSH <- tolower(dis2gene$MeSH_Term)
  dis2gene$MeSH <- trimWhiteSpace(dis2gene$MeSH)
  m.no <- match(gwas.data$disease, dis2gene$disease)
  gwas.data$msh <- dis2gene$MeSH[m.no]
  gwas.data$msh.tree <- dis2gene$MeSH_Tree[m.no]
  gwas.data$msh.cat <- dis2gene$Category[m.no]
  gwas.data$msh.source <- dis2gene$Source[m.no]
  keep.names <- c("snp_id","rsid","PUBMEDID", "pvalue","disease","msh", "msh.tree", "msh.cat")
  rest.names <- names(gwas.data)[!(names(gwas.data) %in% keep.names)]
  gwas.data <- gwas.data[,c(keep.names, rest.names)]
  save(gwas.data, file = "./stopgap_4sources_dbSNP141_msh.RData", compress = TRUE)  
  
  invisible(gwas.data)
}

# Merge GWAS data, LD r2 data and var2gene_vepSimp data together
gwas.ld.var2gene1 <- function(gwas.file="stopgap_4sources_dbSNP141_msh.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData",
                              var2gene.file = "./VarGeneMapping/Var2Gene_vepSimp.RData"){
  
  cat("Start merging GWAS data, LD r2 data and var2gene data together ...", "\n")
  load(gwas.file)
  load(ld.file)
  load(var2gene.file)
  
  names(gwas.data) <- c("snp.gwas.orig", "snp.gwas", "pubmedid","pvalue","disease","msh","msh.id","msh.tree","msh.cat","init.samp", "rep.samp", "source")
  names(ld.snps.r2) <- c("chr.ld", "pos.ld", "ref.ld", "alt.ld", "chr.gwas", "pos.gwas", "snp.gwas","ref.gwas",
                         "alt.gwas", "snp.ld", "r2", "af.1kg","asn.af.1kg", "amr.af.1kg",
                         "afr.af.1kg", "eur.af.1kg")
  ld.snps.r2 <- ld.snps.r2[,c("chr.ld","pos.ld","ref.ld","alt.ld","chr.gwas","pos.gwas","snp.gwas",    
                              "ref.gwas","alt.gwas","snp.ld", "r2", "af.1kg","asn.af.1kg", "amr.af.1kg",
                              "afr.af.1kg", "eur.af.1kg")]
  var2gene.vepSimp <- subset(var2gene.vepSimp, CHR.ld==chr.gencode)
  var2gene.vepSimp$chr.gencode <- NULL
  names(var2gene.vepSimp) <- c("gene","chr.ld","pos.ld", "snp.ld", "cat.rdb", "dhs.type", 
                            "dhs.fpr","chiapet.cell","chiapet.type","fantom5.tissue","fantom5.fpr","fantom5.type","eqtl.ref","eqtl.tissue","eqtl.scoretype","eqtl.score","eqtl.gene.post","eqtl.snp.post","vep.phylop","vep.canonical","vep.conseq","vep.aa","vep.condel","gene.start","gene.end","distance.from.gene")
   
  var2gene.simplified <- var2gene.vepSimp
  
  gwas.r2 <- merge(gwas.data, ld.snps.r2[,c("snp.gwas", "snp.ld","chr.ld", "r2", "af.1kg", "asn.af.1kg", "amr.af.1kg", "afr.af.1kg", "eur.af.1kg")], by.x = "snp.gwas", by.y="snp.gwas", all = TRUE)
 
  system("mkdir GWAS_LD_var2gene")
  save(gwas.r2, file = "./GWAS_LD_var2gene/Mergedata_gwas_r2.RData", compress = TRUE)
 
  chr.no <- c(as.character(1:22), "X")  

#  chr.no <- c(as.character(16:22), "X")
  
  for(i in chr.no){
  
  chr.x1 <- gwas.r2[gwas.r2$chr.ld==i,]
  gwas.r2.var2gene.full <- merge(chr.x1, var2gene.simplified[var2gene.simplified$chr.ld==i,], 
                                 by.x = "snp.ld", by.y="snp.ld", all = TRUE)
 
   gwas.r2.var2gene.withGene <- subset(gwas.r2.var2gene.full, !is.na(gene))
  
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
  #save(gwas.r2.var2gene.full, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_full_r2_0.5.RData", compress = TRUE)
  #save(gwas.r2.var2gene.withGene, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5.RData", compress = TRUE)
  
  gwas.r2.var2gene.full <- subset(gwas.r2.var2gene.full, r2>=0.7)
  save(gwas.r2.var2gene.full, file = paste("./GWAS_LD_var2gene/Mergedata_VEPsimplified_full_r2_0.7_chr",i,".RData",sep=""), compress = TRUE)
  gwas.r2.var2gene.withGene <- subset(gwas.r2.var2gene.withGene, r2>=0.7)
  save(gwas.r2.var2gene.withGene, file = paste("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr",i,".RData",sep=""), compress = TRUE)
  }
  invisible(gwas.r2.var2gene.full)
  
}

# find MHC genes:
find.mhc.genes <- function(gencode.file = "./Data/gencode_v19_clean.RData")
{
  load(gencode.file)
  x <- subset(gencode, chrom %in% "6" & start >= 26016335 & end <= 33548071)
  mhc <- unique(x$Gene)
  
  invisible(mhc)
}

# remove MHC genes from the merged r2=0.7 data:
rm.mhc.genes <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr6.RData",
                         lRmMHC=TRUE){
 
  load(data.file)
  mhc.genes <- find.mhc.genes(gencode.file = "./Data/gencode_v19_clean.RData")
  write.table(mhc.genes, "./GWAS_LD_var2gene/MHC_genes.txt", sep = "\t", row.names = F)
  if (lRmMHC) gwas.r2.var2gene.withGene <- subset(gwas.r2.var2gene.withGene, !(gene %in% trimWhiteSpace(mhc.genes)))  
  
  save(gwas.r2.var2gene.withGene, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr6.RData", compress = TRUE)
  
  invisible(gwas.r2.var2gene.withGene)
}

# Add gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC.RData data
gene.score <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr",
                       severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt"){
  
 cat("Start adding gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7.RData data ...", "\n")
  
 chr.no <- c(as.character(1:22), "X")  
 for(i in chr.no){

  load(paste(data.file,i,".RData",sep=""))
  severity <- read.delim(severity.file,
                         header = TRUE, comment.char = "", sep = "\t", as.is = TRUE)
  gwas.r2.var2gene.withGene$severity <- severity$Severity[match(gwas.r2.var2gene.withGene$vep.conseq, severity$Consequences)]
  
  
  # 1: phylop
  gwas.r2.var2gene.withGene$score1 <- NA
  gwas.r2.var2gene.withGene$score1 <- ifelse(gwas.r2.var2gene.withGene$phylop.score >= 1.5, 1,0)
 
  # 3: cat.rdb
#  cat.lookup <- data.frame(cat=c("1a","1b","1c","1d","1e","1f","2a","2b","2c", "3a","3b", "4"), score=c(rep(1, 9), rep(0, 3)))
#  gwas.r2.var2gene.withGene$score3 <- cat.lookup$score[match(gwas.r2.var2gene.withGene$cat.rdb, cat.lookup$cat)]
  
  # 2: dhs.cor  # Add this one as well
  gwas.r2.var2gene.withGene$dhs.fpr.bin <- as.character(cut(gwas.r2.var2gene.withGene$dhs.fpr, breaks=c(0,0.6,0.85,1), right=FALSE, include.highest=TRUE))
  gwas.r2.var2gene.withGene$dhs.fpr.bin[which((gwas.r2.var2gene.withGene$dhs.fpr) <= 1 & (gwas.r2.var2gene.withGene$dhs.fpr) >= 0.85)] <- "[0.85,1]"
  fpr.lookup <- data.frame(dhs.fpr.bin=c("[0,0.6)", "[0.6,0.85)", "[0.85,1]"), score=c(2,1,0))
  gwas.r2.var2gene.withGene$score2 <- fpr.lookup$score[match(gwas.r2.var2gene.withGene$dhs.fpr.bin, fpr.lookup$dhs.fpr.bin)]

  # 3: deltaSVM
  gwas.r2.var2gene.withGene$deltasvm.bin <- as.character(cut(gwas.r2.var2gene.withGene$deltasvm.quantile.95p, breaks=c(0,5.5,9.7), right=FALSE, include.highest=TRUE))
  gwas.r2.var2gene.withGene$deltasvm.bin[which(gwas.r2.var2gene.withGene$deltasvm.quantile.95p <= max(gwas.r2.var2gene.withGene$deltasvm.quantile.95p,na.rm=T) & gwas.r2.var2gene.withGene$deltasvm.quantile.95p >= 9.7)] <- paste("[9.7,",max(gwas.r2.var2gene.withGene$deltasvm.quantile.95p,na.rm=T),"]",sep="")
  svm.lookup <- data.frame(deltasvm.bin=c("[0,5.5)", "[5.5,9.7)", paste("[9.7,",max(gwas.r2.var2gene.withGene$deltasvm.quantile.95p,na.rm=T),"]",sep="")), score=c(0,1,2))
  gwas.r2.var2gene.withGene$score3 <- svm.lookup$score[match(gwas.r2.var2gene.withGene$deltasvm.bin, svm.lookup$deltasvm.bin)]
  # 4: cato
  gwas.r2.var2gene.withGene$cato.bin <- as.character(cut(gwas.r2.var2gene.withGene$cato.score, breaks=c(0.1,0.25), right=FALSE, include.highest=TRUE))
  gwas.r2.var2gene.withGene$cato.bin[which(gwas.r2.var2gene.withGene$cato.score <= max(gwas.r2.var2gene.withGene$cato.score,na.rm=T) & gwas.r2.var2gene.withGene$cato.score >= 0.25)] <- paste("[0.25,",max(gwas.r2.var2gene.withGene$cato.score,na.rm=T),"]",sep="")
  cato.lookup <- data.frame(cato.bin=c("[0.1,0.25)", paste("[0.25,",max(gwas.r2.var2gene.withGene$cato.score,na.rm=T),"]",sep="")), score=c(1,2))
  gwas.r2.var2gene.withGene$score4 <- cato.lookup$score[match(gwas.r2.var2gene.withGene$cato.bin, cato.lookup$cato.bin)]
  
  # 5: chiapet
  gwas.r2.var2gene.withGene$score5 <- NA
  gwas.r2.var2gene.withGene$score5.1 <- NA

  gwas.r2.var2gene.withGene$score5[!is.na(gwas.r2.var2gene.withGene$chiapet.celltype)] <- 1
  # if presenting in multiple cell lines, further add 1 
  gwas.r2.var2gene.withGene$score5.1[!is.na(gwas.r2.var2gene.withGene$chiapet.celltype) & grepl(";",gwas.r2.var2gene.withGene$chiapet.celltype)] <- 1
  
  # 6: fantom5
  gwas.r2.var2gene.withGene$fantom5.fpr.bin <- as.character(cut(gwas.r2.var2gene.withGene$fantom5.fpr, breaks=c(0,0.6,0.85,1), right=FALSE, include.highest=TRUE))
  gwas.r2.var2gene.withGene$fantom5.fpr.bin[which((gwas.r2.var2gene.withGene$fantom5.fpr) <= 1 & (gwas.r2.var2gene.withGene$fantom5.fpr) >= 0.85)] <- "[0.85,1]"
  fpr.lookup <- data.frame(fantom5.fpr.bin=c("[0,0.6)", "[0.6,0.85)", "[0.85,1]"), score=c(2,1,0))
  gwas.r2.var2gene.withGene$score6 <- fpr.lookup$score[match(gwas.r2.var2gene.withGene$fantom5.fpr.bin, fpr.lookup$fantom5.fpr.bin)]
 
  # 7: eqtl single-tissue
  gwas.r2.var2gene.withGene$score7 <- NA
  gwas.r2.var2gene.withGene$score7[!is.na(gwas.r2.var2gene.withGene$.eqtl.tissue.pval)] <- 2
  gwas.r2.var2gene.withGene$score7.1 <- NA
  gwas.r2.var2gene.withGene$score7.1[!is.na(gwas.r2.var2gene.withGene$single.eqtl.tissue.pval)  & grepl(";",gwas.r2.var2gene.withGene$single.eqtl.tissue.pval)] <- 2
  # 8: eqtl multi-tissue
  gwas.r2.var2gene.withGene$score8 <- NA
  gwas.r2.var2gene.withGene$score8[!is.na(gwas.r2.var2gene.withGene$multi.eqtl.tissue)& is.na(gwas.r2.var2gene.withGene$single.eqtl.tissue.pval)] <- 3
  gwas.r2.var2gene.withGene$score8[!is.na(gwas.r2.var2gene.withGene$single.eqtl.tissue.pval) & !grepl(";",gwas.r2.var2gene.withGene$single.eqtl.tissue.pval) &!is.na(gwas.r2.var2gene.withGene$multi.eqtl.tissue)] <- 2
 
  # 9: vep
  gwas.r2.var2gene.withGene$score9 <- NA
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity <= 5] <- 5
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity <= 10 & gwas.r2.var2gene.withGene$severity >= 6] <- 4
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity <= 17 & gwas.r2.var2gene.withGene$severity >= 11] <- 2
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity <= 33 & gwas.r2.var2gene.withGene$severity >= 18] <- 1
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity == 34] <- 0
  # Severity 6~10: for missense SNPs: +3 if condel deleterious, else neutral +2
  is.neu <- rep(FALSE, nrow(gwas.r2.var2gene.withGene))
  is.neu[grep("neutral", gwas.r2.var2gene.withGene$vep.condel)] <- TRUE
  gwas.r2.var2gene.withGene$score9[gwas.r2.var2gene.withGene$severity <= 10 & gwas.r2.var2gene.withGene$severity >= 6 & gwas.r2.var2gene.withGene$vep.conseq =="missense_variant" & is.neu] <- 3
  
  # 10: CHi-C
  gwas.r2.var2gene.withGene$score10 <- NA
  gwas.r2.var2gene.withGene$score10.1 <- NA

  gwas.r2.var2gene.withGene$score10[!is.na(gwas.r2.var2gene.withGene$chic.tissue)] <- 1
  # if presenting in multiple cell lines, further add 1 
  gwas.r2.var2gene.withGene$score10.1[!is.na(gwas.r2.var2gene.withGene$chic.tissue) & grepl(";",gwas.r2.var2gene.withGene$chic.tissue)] <- 1
 
  
  gwas.r2.var2gene.withGene$var2gene.score <- apply(gwas.r2.var2gene.withGene[,c("score1", "score2", "score3", "score4", "score5","score5.1", "score6", "score7","score7.1","score8","score9","score10","score10.1")], 1, sum, na.rm=TRUE)
  gwas.r2.var2gene.withGene$gene.score <-gwas.r2.var2gene.withGene$r2*gwas.r2.var2gene.withGene$var2gene.score
  
  gwas.r2.var2gene.withGene$dhs.fpr.bin <- NULL
  gwas.r2.var2gene.withGene$deltasvm.bin <- NULL
  gwas.r2.var2gene.withGene$cato.bin <- NULL
  gwas.r2.var2gene.withGene$fantom5.fpr.bin <- NULL
  gwas.r2.var2gene.withGene$score1 <- NULL
  gwas.r2.var2gene.withGene$score2 <- NULL
  gwas.r2.var2gene.withGene$score3 <- NULL
  gwas.r2.var2gene.withGene$score4 <- NULL
  gwas.r2.var2gene.withGene$score5 <- NULL
  gwas.r2.var2gene.withGene$score5.1 <- NULL
  gwas.r2.var2gene.withGene$score6 <- NULL
  gwas.r2.var2gene.withGene$score7 <- NULL
  gwas.r2.var2gene.withGene$score7.1 <- NULL
  gwas.r2.var2gene.withGene$score8 <- NULL
  gwas.r2.var2gene.withGene$score9 <- NULL
  gwas.r2.var2gene.withGene$score10 <- NULL
  gwas.r2.var2gene.withGene$score10.1 <- NULL
  gwas.r2.var2gene.withGene$id <- NULL
  gwas.r2.var2gene.withGene$chiapet.freq <- NULL
  gwas.r2.var2gene.withGene$chr.ld.y <- NULL
  names(gwas.r2.var2gene.withGene)[names(gwas.r2.var2gene.withGene)=="chr.ld.x"] <- "chr.ld"
  # remove repeated snp.gwas+snp.ld+gene names due to # of eqtl tissues
  # aggregate # of eqtl tissues 

  save(gwas.r2.var2gene.withGene, file = paste("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_chr",i,".RData",sep=""), compress = TRUE)
}
 
  invisible(gwas.r2.var2gene.withGene)
  
}
  
#ld.data: Get the LD data as a starting point for SNP locus clustering:
ld.data <- function(data.file="./LocusClustering/gwas_r2.RData"){
  load(data.file)
  system("mkdir LocusClustering")
  gwas.r2.var2gene.full.ld <- gwas.r2[,c("snp.ld","snp.gwas","disease","pvalue","pubmedid","disease","r2","chr.ld", "pos.ld")]
  gwas.r2.var2gene.full.ld <- gwas.r2.var2gene.full.ld[!duplicated(gwas.r2.var2gene.full.ld),]
  gwas.r2.var2gene.full.ld <- gwas.r2.var2gene.full.ld[order(gwas.r2.var2gene.full.ld$disease,
                                                             gwas.r2.var2gene.full.ld$pvalue,
                                                             gwas.r2.var2gene.full.ld$snp.gwas,
                                                             gwas.r2.var2gene.full.ld$snp.ld),]
 chr.no <- c(as.character(1:22), "X")  
 for(i in chr.no){
  ld <- subset(gwas.r2.var2gene.full.ld, chr.ld %in% i)

  #save(gwas.r2.var2gene.full.ld, file = "./LocusClustering/ld.RData", compress = TRUE)
  save(ld, file =paste("./LocusClustering/ld",i,".RData",sep=""), compress = TRUE)
}

  invisible(gwas.r2.var2gene.full.ld)
}

# locus.cluster1: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm does not allow the LD (not GWAS) snps presenting in mulitiple locus clusters
locus.cluster <- function(ld.file="./LocusClustering/ld"){
 
 chr.no <- c(as.character(1:22), "X")  
 for(k in chr.no){

  load(paste(ld.file,k,".RData",sep=""))
  
  gwas.r2.var2gene.full.r2 <- ld
  
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
  ld.clusters <- do.call("rbind", clusters)
  save(ld.clusters, file = paste("./LocusClustering/STOPGAP_SNP_Cluster_chr",k,".RData",sep=""), compress = TRUE)
  }
}

#Fine mapping of casual variants using PICS 
fine.map.pics <-function(cluster.file = "./LocusClustering/STOPGAP_SNP_Cluster_",
                         gwas.file="stopgap_4sources_dbSNP141_msh.RData"){

  cat("Start merging cluster data and the merged data before bestLD data together ...", "\n")

  chr.no <- c(as.character(1:22),"X")  

  load(gwas.file)
  gwas.data$id=paste(gwas.data$disease,gwas.data$snp_id,gwas.data$PUBMEDID,gwas.data$pvalue,sep="_")
  gwas.data = gwas.data[,c("id","Source")]
  source.yes <- c("PheWAS","nhgri")

  for(i in chr.no){

  load(paste(cluster.file,i,".RData",sep=""))
  stopgap.snp.cluster$id=paste(stopgap.snp.cluster$disease,stopgap.snp.cluster$snp.gwas,stopgap.snp.cluster$pubmedid,stopgap.snp.cluster$pvalue,sep="_")
#  names(stopgap.snp.cluster)[names(stopgap.snp.cluster)=="start.gwassnp"] <- "clust.init.gwassnp"
  stopgap.snp.cluster <-   merge(stopgap.snp.cluster, gwas.data,by ="id")
  stopgap.snp.cluster$id=NULL
  source.list <- names(table(stopgap.snp.cluster$Source))
  source.no <- source.list[!source.list %in% source.yes]
  source.yes <- source.list[source.list %in% source.yes]

 all= subset(stopgap.snp.cluster, Source %in% source.yes)
 for(j in source.no){ 
  ngc <- subset(stopgap.snp.cluster, Source %in% j)
 # fine mapping step
  ngc <- ngc[order(ngc$pubmedid,
                 ngc$disease,
                 ngc$cluster,
                 ngc$pvalue),]
  ngc$id <- paste(ngc$pubmedid, ngc$disease, ngc$cluster, sep="_")
  ngc$id1 <- paste(ngc$id, ngc$pvalue, sep="_")
  
  # Add the p-value rank
  tmp <- ngc[,c("id", "id1")]
  tmp1 <- tmp[match(unique(tmp$id1),tmp$id1),]
  tmp2 <- tmp1[match(unique(tmp1$id),tmp1$id),]

  id.count <- count(tmp1, "id")
  # Put back the order. 
  id.count <- id.count[match(tmp2$id, id.count$id),]
  tmp1$rank <- unlist(lapply(id.count$freq, function(x) 1:x))
  ngc$gwassnp.p.rank <- tmp1$rank[match(ngc$id1, tmp1$id1)]
  ngc <- subset(ngc, gwassnp.p.rank == 1)
  ngc$id <- NULL
  ngc$id1 <- NULL
  ngc$gwassnp.p.rank <- NULL
  rm(tmp)
  rm(tmp1)
  all=rbind(all,ngc)
}
  stopgap.snp.cluster=all
  stopgap.snp.cluster$pics=NULL
  stopgap.snp.cluster$Source=NULL
  stopgap.snp.cluster$id=NULL

  # calculate PICS
   stopgap.snp.cluster$pvalue=as.double(  stopgap.snp.cluster$pvalue)
    stopgap.snp.cluster$pvalue=ifelse(  stopgap.snp.cluster$pvalue == 0, 1e-250,   stopgap.snp.cluster$pvalue)
  
    stopgap.snp.cluster$SD=ifelse(!is.na(  stopgap.snp.cluster$r2),sqrt(1-sqrt(  stopgap.snp.cluster$r2)^6.4)*sqrt(-log10(  stopgap.snp.cluster$pvalue))/2,0)
    stopgap.snp.cluster$Mean=ifelse(!is.na(  stopgap.snp.cluster$r2),  stopgap.snp.cluster$r2*-log10(  stopgap.snp.cluster$pvalue),1-log10(  stopgap.snp.cluster$pvalue))
  
  #  stopgap.snp.cluster$Stat=(-log10(  stopgap.snp.cluster$pvalue)-  stopgap.snp.cluster$Mean)/  stopgap.snp.cluster$SD
     stopgap.snp.cluster$prob=ifelse(  stopgap.snp.cluster$SD==0,1,1-pnorm(-log10(  stopgap.snp.cluster$pvalue),  stopgap.snp.cluster$Mean,  stopgap.snp.cluster$SD))  
  # normalized the probabilities so that the total of their probability summed to 1 
     stopgap.snp.cluster$id=paste(  stopgap.snp.cluster$cluster,  stopgap.snp.cluster$disease,  stopgap.snp.cluster$pubmedid,sep="_")
   prob.sum=aggregate(prob~id,  stopgap.snp.cluster,sum)
   colnames(prob.sum)=c("id","prob_sum")
     stopgap.snp.cluster=merge(  stopgap.snp.cluster,prob.sum,by="id")
     stopgap.snp.cluster$pics=  stopgap.snp.cluster$prob/  stopgap.snp.cluster$prob_sum
   # Clean up some columns
   stopgap.snp.cluster=stopgap.snp.cluster[,-match(c("SD","Mean","prob","prob_sum"), names(stopgap.snp.cluster))]

  stopgap.snp.cluster <- droplevels(all)

  save(stopgap.snp.cluster, file = paste( "./LocusClustering/STOPGAP_SNP_Cluster_chr",i,".RData",sep=""), compress = TRUE)
 }
}


# Create a function genomics and annotation evidence column 
sg.evid <- function(data) {
  estr <- rep("", nrow(data))
 
  lVep <- !is.na(data$vep.severity) & data$vep.severity <=33
  lNS <- !is.na(data$vep.condel)
  x <- substr(as.character(data$vep.condel[lNS]), 1, 1)
  estr[lVep] <- paste(estr[lVep], rep("V", sum(lVep)), sep = "")
  estr[lVep & lNS] <- paste(estr[lVep & lNS], x, sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  leQTL <- !is.na(data$single.eqtl.tissue.pvalue) | !is.na(data$multi.eqtl.tissue)
  estr[leQTL] <- paste(estr[leQTL], rep("E", sum(leQTL)), sep = "")
 
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  lChiC <- !is.na(data$chic.tissue)
  estr[lChiC] <- paste(estr[lChiC], rep("CC", sum(lChiC)),
                      data$chic.type[lChiC],
                      sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  lChia <- !is.na(data$chiapet.celltype)
  estr[lChia] <- paste(estr[lChia], rep("CP", sum(lChia)),
                       data$chiapet.type[lChia],
                       sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  lDhs <- !is.na(data$dhs.fpr) & data$dhs.fpr < 0.85
  estr[lDhs] <- paste(estr[lDhs], rep("D", sum(lDhs)), data$dhs.type[lDhs],
                      sep = "")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  lF5 <- !is.na(data$fantom5.fpr) & data$fantom5.fpr < 0.85
  estr[lF5] <- paste(estr[lF5], rep("F5", sum(lF5)), sep = "")

  estr <- paste(estr, rep(";", length(estr)), sep = "")

  lSvm <- !is.na(data$deltasvm.quantile.95p) & data$deltasvm.quantile.95p >= 5.5
  estr[lSvm] <- paste(estr[lSvm], rep("DS", sum(lSvm)), sep="")
  estr <- paste(estr, rep(";", length(estr)), sep = "")
  
  lCat <- !is.na(data$cato.score) & data$cato.score >= 0.1
  estr[lCat] <- paste(estr[lCat], rep("CT", sum(lCat)), sep="")
  estr <- paste(estr, rep(";", length(estr)), sep = "")
 
  lPhy <- !is.na(data$phylop.score) & data$phylop.score >= 1.5
  estr[lPhy] <- paste(estr[lPhy], rep("P", sum(lPhy)), sep="")
  
  estr <- paste(estr, rep(";", length(estr)), sep = "")

  
  return(factor(estr))
}

#This is grouped by publication-trait-GWAS SNP-gene, i.e. we reduce each reported genetic association to the best variant-to-gene hypothesis for each possible gene mapping.
# For each unique GWAS SNP, rank its LD SNPs' genes based on their gene scores and r2 if gene scores tie
gene.rank <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_chr"){
  
  cat("For each unique GWAS SNP, start ranking its LD SNPs' genes based on their gene scores...", "\n")

  chr.no <- c(as.character(1:22), "X")  
  for(i in chr.no){

  load(paste(data.file,i,".RData",sep=""))
  
  data.clean <-gwas.r2.var2gene.withGene
  data.clean <- subset(data.clean, !is.na(pvalue))
  data.clean$vep.canonical[which(data.clean$vep.canonical=="-")] <- NA
  data.clean$vep.aa[which(data.clean$vep.aa=="-")] <- NA
  data.clean$vep.condel[which(data.clean$vep.condel=="-")] <- NA
  names(data.clean)[names(data.clean)=="severity"] <- "vep.severity"

  data.clean <- data.clean[!is.na(data.clean$gene.score),]
   
  data.clean <- data.clean[order(data.clean$snp.gwas,
                                -data.clean$gene.score),]
  data.clean$id <- paste(data.clean$snp.gwas, data.clean$gene, sep="_")
  #data.clean.a <- subset(data.clean, snp.gwas %in% c("rs13111494","rs884304", "rs4665058"))
  data.a <- data.clean[, c("snp.gwas", "gene")]
  data.a1 <- data.a[!duplicated(data.a),]
  data.a1$id <- paste(data.a1$snp.gwas, data.a1$gene, sep="_")
  
  library(plyr)
  count1 <- count(data.a, c("snp.gwas", "gene"))  
#  snp.gwas       gene freq
# rs9982597 AP000295.9    2 # of pairs
# rs9982597       HUNK    2
# rs9982597     IL10RB    2
  # Please note count() automatically sorts the snp.gwas-gene pair names, so need to put the sorted order back
  count2 <- count(count1[,-3], "snp.gwas")    # Need to remove the freq column (-3), otherwise the numbers are not what I want
#  snp.gwas  freq  
# rs9982597	    3 # of genes 
  count1$id <- paste(count1$snp.gwas, count1$gene,sep="_")
  count1 <- count1[match(data.a1$id, count1$id),]
  count1$rank <- unlist(lapply(count2$freq, function(x) 1:x))
  #data.a$gene <- as.character(data.a$gene)
  #data$snp.gwas <- as.character(data$snp.gwas)
  
  data.clean$gene.rank.random <- count1$rank[match(data.clean$id, count1$id)]
  
  # Rank based on rank(x, ties.method = "max") and rank(x, ties.method = "min")
  data1 <- data.clean[,c("snp.gwas", "gene", "id", "gene.score")]
  data2 <- data1[match(count1$id, data1$id),] 
  data2$genescorer2 <- -1 * data2$gene.score
  # rank 1 means the best gene score and r2
  gene.rank1 <- aggregate(genescorer2 ~ snp.gwas, data2, function(x){rank(x, ties.method = "max")})
  # Divide out ranks with more than one rank (replace by their max in a set of ties
  rank.max <- as.numeric(unlist(gene.rank1$genescorer2))
  count1$rank.max <- rank.max  
  data.clean$gene.rank.max <- count1$rank.max[match(data.clean$id, count1$id)]
  
  gene.rank2 <- aggregate(genescorer2 ~ snp.gwas, data2, function(x){rank(x, ties.method = "min")})
  # Divide out ranks with more than one rank (replace by their min
  rank.min <- as.numeric(unlist(gene.rank2$genescorer2))
  count1$rank.min <- rank.min  
  data.clean$gene.rank.min <- count1$rank.min[match(data.clean$id, count1$id)]  
  
  # Remove the gene.rank.random column
  data.clean$gene.rank.random <- NULL
  data.clean$id <- NULL
  stopgap.data <- data.clean
  
  # add the evidence ***
  stopgap.data$evidence <- sg.evid(stopgap.data)
 
  # add the best gene and evidence given snp.gwas
  best.gene.score <- stopgap.data[,c("snp.gwas", "gene.score", "gene", "evidence")]
  best.gene.score <- best.gene.score[order(best.gene.score$snp.gwas, 
                                               -best.gene.score$gene.score),]
  best.gene.score <- best.gene.score[match(unique(best.gene.score$snp.gwas), best.gene.score$snp.gwas),]
  best.gene.score1 <- best.gene.score[,c("snp.gwas","gene")]
  stopgap.data$best.gene <- best.gene.score1$gene[match(stopgap.data$snp.gwas, best.gene.score1$snp.gwas)]  
  best.gene.score1 <- best.gene.score[,c("snp.gwas","evidence")]
  stopgap.data$best.evidence <- best.gene.score1$evidence[match(stopgap.data$snp.gwas, best.gene.score1$snp.gwas)]  

  save(stopgap.data, file = paste("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_chr",i,".RData",sep=""), compress = TRUE)
  
  # Select the best ld SNPs per gene and each original GWAS SNP based their gene score and r2
  # Order by the snp.gwas, rev(gene.score) and rev(r2), then pick the first row
  
  data.bestld <- stopgap.data
  data.bestld$id <- paste(data.bestld$snp.gwas, data.bestld$pubmedid, data.bestld$disease, data.bestld$gene, sep="_")
  data.bestld <- data.bestld[match(unique(data.bestld$id), data.bestld$id),]
  data.bestld$id <- NULL

  
  stopgap.data.bestld <- data.bestld
  save(stopgap.data.bestld, file = paste("./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_cluster_generank_bestLD_chr",i,".RData",sep=""), compress = TRUE)

 }  
}  

# Add cluser number data to the best LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_chr",
                        cluster.file = "./LocusClustering/STOPGAP_SNP_Cluster_"){
  
  cat("Start merging cluster data and the merged data before bestLD data together ...", "\n")

  chr.no <- c(as.character(1:22), "X")  
  for(i in chr.no){

  load(paste(data.file,i,".RData",sep=""))
  load(paste(cluster.file,i,".RData",sep=""))
  names(stopgap.snp.cluster)[names(stopgap.snp.cluster)=="start.gwassnp"] <- "clust.init.gwassnp"
  stopgap.data.bestld <-   merge(stopgap.data.bestld, stopgap.snp.cluster,
                          by =c("snp.ld","snp.gwas","chr.ld","pos.ld","pvalue","r2","disease","pubmedid"), all.x = TRUE, all.y=FALSE)

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
  library(plyr)
  id.count <- count(tmp1, "id")
  # Put back the order. 
  id.count <- id.count[match(tmp2$id, id.count$id),]
  tmp1$rank <- unlist(lapply(id.count$freq, function(x) 1:x))
  stopgap.data.bestld$gwassnp.p.rank <- tmp1$rank[match(stopgap.data.bestld$id1, tmp1$id1)]
  stopgap.data.bestld$id <- NULL
  stopgap.data.bestld$id1 <- NULL
  rm(tmp)
  rm(tmp1)

  save(stopgap.data.bestld, file = paste( "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster_chr",i,".RData",sep=""), compress = TRUE)
}
  
}

#combine stopgap.data.bestld datasets
Combine.bestld.data <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster_chr"){
  
 cat("Start subset the bestLD data per each gene and msh combination ...", "\n")
 chr.no <- c(as.character(1:22),"X")  
 
 all.result=list()
 k=0
 for(i in chr.no) {
 k<-k+1
  load(paste(data.file,i,".RData",sep=""))
  all.result[[k]]=stopgap.data.bestld
 }
 stopgap.data.bestld=do.call("rbind", all.result)
 save(stopgap.data.bestld, file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData", compress = TRUE)

}


# Generate the bestld data per gene-MeSH.
# Rank by the p-value bins and then by gene scores for  bestld data per gene-MeSH.
# P-value bins (-log10(P)): > 12 as 1; 8-12 as 2 and <8 as 3.
mesh.gene <- function(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster.RData"){
  
 cat("Start subset the bestLD data per each gene and msh combination ...", "\n")
 
 load(data.file)
  
  stopgap.data.bestld <- subset(stopgap.data.bestld, !is.na(msh))
  # ~ 0.25% rows get removed due to missing msh terms
  
  stopgap.data.bestld$id <- paste(stopgap.data.bestld$gene, stopgap.data.bestld$msh, sep="_")
  
  # p-value bins
  stopgap.data.bestld$p.cat <- cut(as.double(stopgap.data.bestld$pvalue), breaks=c(0, 1e-12, 1e-8, 1)) 
  lookup <- data.frame(cat=c("(0,1e-12]", "(1e-12,1e-08]", "(1e-08,1]"), bin=c(1,2,3))
  stopgap.data.bestld$p.bin <- lookup$bin[match(stopgap.data.bestld$p.cat, lookup$cat)]
  
  # Rank the gene_msh by first p-value and then gene score (decreasing)
  stopgap.data.bestld <- stopgap.data.bestld[order(stopgap.data.bestld$id,
                                                   stopgap.data.bestld$p.bin,
                                                   -stopgap.data.bestld$gene.score),]
  #subset(stopgap.data.bestld, id=="A1CF_hyperuricemia")
  library(plyr)
  
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
  
  save(stopgap.bestld.gene.mesh, file = "./Mesh_gene/STOPGAP_r2_0.7_MHC_bestLD_gene_mesh.RData", compress = TRUE)
#}  
  invisible(stopgap.bestld.gene.mesh)
  
}

# merge.orphanet: Merge stopgap.gene.mesh with latest Orphanet data processed on 11/13/2104
# Create merge between Orphanet, and stopgap.gene.mesh
# Change "Link" to "pubmedid" and add the OrphID as "pubmedid" when the Source=="Orphanet". 
# Set Rank=1, pval2 as 0 and GeneScore=max(stopgap.gene.mesh$gene.score)
# Remove MSH.TOP and OrphID columms and change the names to being consistent with stopgap.gene.mesh
merge.orphanet <- function(gm.file = "./Mesh_gene/STOPGAP_r2_0.7_rmMHC_bestLD_gene_mesh.RData",
                       orphanet.file = "./Data/orphanet.RData",
					   gencode.file="./Data/gencode_v19_clean.RData"
					   ){  
  load(gm.file)
  load(orphanet.file)
  load(gencode.file)
  
  
  new.gene.name <- read.delim("./Data/hg19tohg38.txt",header=T)
  names(new.gene.name)<-c("gene.v19","gene")
  new.gene.name<-new.gene.name[!duplicated(new.gene.name),]
 
  stopgap.gene.mesh <- stopgap.bestld.gene.mesh
  orphanet$var2gene.score <- max(stopgap.bestld.gene.mesh$var2gene.score)
  orphanet$gene.score <- max(stopgap.bestld.gene.mesh$gene.score)
  orphanet$gene.rank.max <- 1
  orphanet$gene.rank.min <- 1
  orphanet$pvalue <- 0
  orphanet$best.gene <- omim$gene

  orphanet <- orphanet[,c("disease","pubmedid","gene",         
  "msh","msh.tree","msh.cat","var2gene.score",       
  "gene.score","best.gene","gene.rank.max","gene.rank.min","pvalue","source")]
  
  stopgap.gene.mesh$pubmedid=as.character(stopgap.gene.mesh$pubmedid)
  orphanet$pubmedid=as.character(orphanet$pubmedid)
 
  stopgap.gene.mesh <- sbind(stopgap.gene.mesh, orphanet)
  
  stopgap.gene.mesh <- rename(stopgap.gene.mesh, c(best.gene= 'gene.best', best.evidence="evidence.best",
                                                   best.gene.score='gene.score.best', distance.from.gene="gene.distance",
					                               chiapet.celltype = 'chiapet.cell', deltasvm.top5.celltype="deltasvm.cell",
												   cato.celltype="cato.cell", var2gene.score="v2g.score"))
  stopgap.gene.mesh$eqtl.tissue.pvalue=ifelse((!is.na(stopgap.gene.mesh$single.eqtl.tissue.pvalue) & !is.na(stopgap.gene.mesh$multi.eqtl.tissue)), paste(stopgap.gene.mesh$single.eqtl.tissue.pvalue,"multi-tissues",sep=";"),ifelse((!is.na(stopgap.gene.mesh$single.eqtl.tissue.pvalue) & is.na(stopgap.gene.mesh$multi.eqtl.tissue)),stopgap.gene.mesh$single.eqtl.tissue.pvalue, ifelse((is.na(stopgap.gene.mesh$single.eqtl.tissue.pvalue) & !is.na(stopgap.gene.mesh$multi.eqtl.tissue)),"multi-tissues",NA)))
  stopgap.gene.mesh$single.eqtl.gene.type=NULL
  stopgap.gene.mesh$multi.eqtl.tissue=NULL
  stopgap.gene.mesh$single.eqtl.tissue.pvalue=NULL
  stopgap.gene.mesh$vep.canonical=NULL
  stopgap.gene.mesh$chiapet.type=NULL
  stopgap.gene.mesh$chic.type=NULL
  stopgap.gene.mesh$post=NULL
#  stopgap.gene.mesh$dhs.type=NULL
  stopgap.gene.mesh$max.gene.score=NULL
  stopgap.gene.mesh$fantom5.type=NULL
 # stopgap.cat.rdb=NULL
  
  stopgap.gene.mesh$msh.cat=ifelse(trimWhiteSpace(stopgap.gene.mesh$msh.cat)=="Digestive System","Digestive system", trimWhiteSpace(stopgap.gene.mesh$msh.cat))

  stopgap.gene.mesh$chr.ld=as.factor(stopgap.gene.mesh$chr.ld)
  stopgap.gene.mesh$pvalue=as.double(stopgap.gene.mesh$pvalue)
  stopgap.gene.mesh$af.1kg = as.numeric(stopgap.gene.mesh$af.1kg)
  stopgap.gene.mesh$asn.af.1kg=as.numeric(stopgap.gene.mesh$asn.af.1kg)
  stopgap.gene.mesh$amr.af.1kg=as.numeric(stopgap.gene.mesh$amr.af.1kg)
  stopgap.gene.mesh$afr.af.1kg=as.numeric(stopgap.gene.mesh$afr.af.1kg)

  stopgap.gene.mesh <- rename(stopgap.gene.mesh, c(gene="gene.v19"))
  stopgap.gene.mesh <- merge(stopgap.gene.mesh,new.gene.name,by="gene.v19",all.x=T,all.y=F)
  
  stopgap.gene.mesh$gene <- as.character(stopgap.gene.mesh$gene)
  stopgap.gene.mesh$gene <- ifelse(!is.na(stopgap.gene.mesh$gene), stopgap.gene.mesh$gene, stopgap.gene.mesh$gene.v19)
  stopgap.gene.mesh <- rename(stopgap.gene.mesh, c(gene.best="gene.best.v19"))

  new.gene.name <- rename(new.gene.name, c(gene="gene.best", gene.v19="gene.best.v19"))

  stopgap.gene.mesh <- merge(stopgap.gene.mesh,new.gene.name,by="gene.best.v19",all.x=T,all.y=F)
  
  
  stopgap.gene.mesh$gene.best <- as.character(stopgap.gene.mesh$gene.best)
  stopgap.gene.mesh$gene.best <- ifelse(!is.na(stopgap.gene.mesh$gene.best), stopgap.gene.mesh$gene.best, stopgap.gene.mesh$gene.best.v19)
 
  
var=c("gene","gene.v19","disease","msh","pvalue","pubmedid","source","snp.gwas","snp.ld","r2","evidence","v2g.score",
"gene.score","gene.rank.max","gene.rank.min","gene.distance","gene.best","evidence.best","gene.score.best",
"cluster","vep.conseq","vep.severity","vep.aa","vep.condel","eqtl.tissue.pvalue","chic.tissue",
"chiapet.cell","dhs.type","dhs.fpr","fantom5.tissue","fantom5.fpr","deltasvm.95p","deltasvm.cell","cato.score","cato.motifname",
"cato.cell","vep.phylop","chr.ld","pos.ld","snp.gwas.orig","init.samp","rep.samp","msh.tree",
"msh.cat","af.1kg","asn.af.1kg","amr.af.1kg","afr.af.1kg","eur.af.1kg","gene.start","gene.end",
"clust.init.gwassnp","pics","gwassnp.p.rank", 
"asso.count","num.cluster","num.gwas.snp",
"num.ld.snp","num.gene.cluster")

  stopgap.gene.mesh <- stopgap.gene.mesh[,var]

 stopgap.gene.mesh$pubmedid=ifelse(stopgap.gene.mesh$source == "PheWAS","20335276",stopgap.gene.mesh$pubmedid)
 write.table(stopgap.gene.mesh,
              "stopgap.gene.mesh.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")  
  save(stopgap.gene.mesh, file = "stopgap.gene.mesh.RData", compress = TRUE)   
  
  invisible(stopgap.gene.mesh)
  
}

revice.bestld <- function(gm.file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster.RData",
                          gencode.file="./Data/gencode_v19_clean.RData" ){  
  load(gm.file)
  load(gencode.file)

  new.gene.name <- read.delim("./Data/hg19tohg38.txt",header=T)
  names(new.gene.name)<-c("gene.v19","gene")
  new.gene.name<-new.gene.name[!duplicated(new.gene.name),]
  
 stopgap.bestld <- stopgap.data.bestld
 
 stopgap.bestld <- rename(stopgap.bestld, c(best.gene= 'gene.best', best.evidence="evidence.best",
                                                   best.gene.score='gene.score.best', distance.from.gene="gene.distance",
					                               chiapet.celltype = 'chiapet.cell', deltasvm.top5.celltype="deltasvm.cell",
												   cato.celltype="cato.cell", var2gene.score="v2g.score"))
  stopgap.bestld$eqtl.tissue.pvalue=ifelse((!is.na(stopgap.bestld$single.eqtl.tissue.pvalue) & !is.na(stopgap.bestld$multi.eqtl.tissue)), paste(stopgap.bestld$single.eqtl.tissue.pvalue,"multi-tissues",sep=";"),ifelse((!is.na(stopgap.bestld$single.eqtl.tissue.pvalue) & is.na(stopgap.bestld$multi.eqtl.tissue)),stopgap.bestld$single.eqtl.tissue.pvalue, ifelse((is.na(stopgap.bestld$single.eqtl.tissue.pvalue) & !is.na(stopgap.bestld$multi.eqtl.tissue)),"multi-tissues",NA)))
  stopgap.bestld$single.eqtl.gene.type=NULL
  stopgap.bestld$multi.eqtl.tissue=NULL
  stopgap.bestld$single.eqtl.tissue.pvalue=NULL
  stopgap.bestld$vep.canonical=NULL
  stopgap.bestld$chiapet.type=NULL
  stopgap.bestld$chic.type=NULL
  stopgap.bestld$post=NULL
#  stopgap.bestld$dhs.type=NULL
  stopgap.bestld$max.gene.score=NULL
  stopgap.bestld$fantom5.type=NULL
 # stopgap.cat.rdb=NULL
  
  stopgap.bestld$msh.cat=ifelse(trimWhiteSpace(stopgap.bestld$msh.cat)=="Digestive System","Digestive system", trimWhiteSpace(stopgap.bestld$msh.cat))

  stopgap.bestld$chr.ld=as.factor(stopgap.bestld$chr.ld)
  stopgap.bestld$pvalue=as.double(stopgap.bestld$pvalue)
  stopgap.bestld$af.1kg = as.numeric(stopgap.bestld$af.1kg)
  stopgap.bestld$asn.af.1kg=as.numeric(stopgap.bestld$asn.af.1kg)
  stopgap.bestld$amr.af.1kg=as.numeric(stopgap.bestld$amr.af.1kg)
  stopgap.bestld$afr.af.1kg=as.numeric(stopgap.bestld$afr.af.1kg)

  stopgap.bestld <- rename(stopgap.bestld, c(gene="gene.v19"))
  stopgap.bestld <- merge(stopgap.bestld,new.gene.name,by="gene.v19",all.x=T,all.y=F)
  stopgap.bestld$gene <- as.character(stopgap.bestld$gene)
  stopgap.bestld$gene <- ifelse(!is.na(stopgap.bestld$gene), stopgap.bestld$gene, stopgap.bestld$gene.v19)
  stopgap.bestld <- rename(stopgap.bestld, c(gene.best="gene.best.v19"))

  new.gene.name <- rename(new.gene.name, c(gene="gene.best", gene.v19="gene.best.v19"))

  stopgap.bestld <- merge(stopgap.bestld,new.gene.name,by="gene.best.v19",all.x=T,all.y=F)
  stopgap.bestld$gene.best <- as.character(stopgap.bestld$gene.best)
  stopgap.bestld$gene.best <- ifelse(!is.na(stopgap.bestld$gene.best), stopgap.bestld$gene.best, stopgap.bestld$gene.best.v19)
 
 
var=c("gene","gene.v19","disease","msh","pvalue","pubmedid","source","snp.gwas","snp.ld","r2","evidence","v2g.score",
"gene.score","gene.rank.max","gene.rank.min","gene.distance","gene.best","evidence.best","gene.score.best",
"cluster","vep.conseq","vep.severity","vep.aa","vep.condel","eqtl.tissue.pvalue","chic.tissue",
"chiapet.cell","dhs.type","dhs.fpr","fantom5.tissue","fantom5.fpr","deltasvm.95p","deltasvm.cell","cato.score","cato.motifname",
"cato.cell","vep.phylop","chr.ld","pos.ld","snp.gwas.orig","init.samp","rep.samp","msh.tree",
"msh.cat","af.1kg","asn.af.1kg","amr.af.1kg","afr.af.1kg","eur.af.1kg","gene.start","gene.end",
"clust.init.gwassnp","pics","gwassnp.p.rank")

  stopgap.bestld <- stopgap.bestld[,var]
  write.table(stopgap.bestld,
              "stopgap.bestld.txt",
              sep = "\t", row.names = FALSE, quote = F, na = "")  
  save(stopgap.bestld, file = "stopgap.bestld.RData", compress = TRUE)   
  
  invisible(stopgap.bestld)
  
}
