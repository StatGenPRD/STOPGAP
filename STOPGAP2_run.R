#--------Run STOPGAP2 R functions to generate STOPGAP2 datasets -----------------#
# Author: Judong Shen and Kijoung Song

setwd("STOPGAP.working.directory")
library(gplots)
library(plyr)
library(reshape)
#library(cba)
#library(ggplot2)

source("STOPGAP2_functions.R")
# (1) merge the grasp, nhgri and gwasdb datasets
p.threshold <- 1E-04
grasp.import(file = "./Data/GraspFullDataset2.zip", p.thres=p.threshold)
NHGRI.import(file = "./Data/gwas_catalog_v1.0-associations_e85_r2016-08-01.tsv", p.thres=p.threshold)
gwasdb.import(file = "./Data/gwasdb_20150819_snp_trait.gz", p.thres=p.threshold)
phewas.import(file = "./Data/phewas-catalog.csv", p.thres=1E-04)
data.merge(nhgri, grasp, gwasdb, phewas)

# Identify the unique GWAS SNPs 
gwas.snps(gwas.file="stopgap_4sources_clean.RData")

# Coordinate lookup and LD calculation based on the 1KG data for the new SNPs identified in this version
# if you want to calcuate LD SNPs based on your data, please see 'README.txt'.
#run.ld()

# Update rsIDs in the gwas data (data) to dbSNP141 version
rsID.update(gwas.file="stopgap_4sources_clean.RData", 
                        coor.file="./STOPGAP2_LDResults/rsID_Coordinates.txt")


# Find all variants in LD with GWAS variants at r2 > 0.5 
ld.snps(gwas.data="stopgap_4sources_dbSNP141.RData",
                    ld.path="./STOPGAP2_LDResults")
  
# updata.ld.data: updata ld data:
# Remove SNP.ld with two positions in both ld.snps and ld.snps.r2
# Use chr:pos to replace the missing ld SNP ids in both ld.snps and ld.snps.r2
updata.ld.data(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.RData", 
                ld.r2.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.RData")


# Update the ld.snps.r2 and ld_snps data by including the GWAS SNPs data which don't have any LD information there
update.ld.snps.r2(gwas.file="stopgap_4sources_dbSNP141.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.RData",
                              rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt")

# Update the ld.snps.r2 and ld_snps data by removing the LD SNPs more than 500kb away from the GWAS SNPs,
# If any GWAS SNPs get deleted, put them back by including r2=1 there
update1.ld.snps.r2(gwas.file="stopgap_4sources_dbSNP141.RData",
                               ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.RData",
                               dis = 1000)

# Update the ld.snps.r2 data by adding the AF information from 1KG phase I data ...
r2.frq(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.RData",
                   frq.file="./Data/1KG_AF/gws.frq.RData")

# For each LD SNP, (1) identify each dhscor where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.dhscor(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                      dhscor.file="./Data/genomic_resource/dhscor.RData",
                      mpos.rng.p = c("Chr.pdhs", "Start.pdhs", "Stop.pdhs"),
                      mpos.rng.d = c("Chr.ddhs", "Start.ddhs", "Stop.ddhs"))

# ----- Merge the CHIA-PET data
# ld.chiapet: For each LD SNP, (1) indentify each chiapet where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.chiapet(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                       chiapet.file="./Data/genomic_resource/CHIA.PET_5celltypes.RData",
                       mpos.rng.d = c("dChr", "dStart", "dEnd"))

# ld.eQTL: For each LD SNP, (1) indentify each eQTL where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.eQTL(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                    eQTL.file="./Data/genomic_resource/eQTL_single_44Tissues.RData")

ld.multi.eQTL(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData" ,
                          eQTL.multi.file="./Data/genomic_resource/eQTL_multi_9Tissues.RData")


# ld.fantom5: For each LD SNP, (1) indentify each fantom5 where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.fantom5(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",  
                       fantom5.file="./Data/genomic_resource/FANTOM5.RData",
                       mpos.rng.d = c("chr", "startpos", "endpos"))

# ld.CHIC: For each LD SNP, (1) indentify each CHIC where SNP falls within promoter or distal genomic region; 
#                           (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.CHIC(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                       CHIC.file="./Data/genomic_resource/chic.RData",
                       mpos.rng.d = c("chr", "startpos", "endpos"))

# add cato to ld snps 
ld.cato(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                   cato.file="./Data/genomic_resource/cato.RData")

# gwas.ld.deltaSVM.RData is generated and downloaded at Data/genomic_resource 
ld.deltaSVM()

# gwas.ld.phylop.RData is generated and downloaded at Data/genomic_resource 
ld.phylop()

# ld.Posgene: For each LD SNP, (1) identify genes falling in the 5Kb windows (upstream 5kb, downstream 5kb) 
ld.posgene(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                       Posgene.file="/Data/genomic_resource/gencode_v19_clean.RData")

# ld.vep: Merge gwas ld SNPs and VEP information. 
# gwas.ld.vep.RData is generated and downloaded at Data/genomic_resource 
ld.vep()


# ld.vep.simplify: Simplify the full version of merged gwas ld SNPs and VEP data. 
ld.vep.simplify(vep.file="./VarGeneMapping/gwas.ld.vep.full.RData", 
                           severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")

#merge.ld.simp: Merge based on SNP.ld and gene from each data
# % of unique (ld) GWAS SNPs with at least one genes
merge.ld.simp(dhscor.file ="./VarGeneMapping/gwas.ld.dhscor.RData",
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

# Subset of the var2gene data by using the union of RefSeq and Gencode gene names
# Also Update the var2gene mapping data by removing any rows with conflicting chr.ld with those from the chr information from gencode data
var2gene.filter(var2gene.file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData",
                            gencode.file = "./Data/gencode_v19_clean.RData")

# gwas.mesh: add mesh terms to the GWAS data. 
# Run it on Windows only
gwas.mesh(gwas.file="./stopgap_4sources_dbSNP141.RData",
                      disease2mesh.file = "./Data/MeSH_Aug2016.txt")

# Merge GWAS data, LD r2 data and var2gene_vepSimp data together
gwas.ld.var2gene1(gwas.file="stopgap_4sources_dbSNP141_msh.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData",
                              var2gene.file = "./VarGeneMapping/Var2Gene_vepSimp.RData")

rm.mhc.genes(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr6.RData",
                         lRmMHC=TRUE)

# Add gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7.RData data
gene.score(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_chr",
                       severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")

# gene.score1(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5_rmMHC.RData",
#             severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")
# 
# gene.rank1(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5_rmMHC_genescore.RData")

#ld.data: Get the LD data as a starting point for SNP locus clustering:
ld.data(data.file="./LocusClustering/gwas_r2.RData")

# locus.cluster1: clustering all the GWAS SNPs and LD SNPs into locus clusters based on LD r2 information between the GWAS SNPs and LD SNPs
# Note the ld data in the ld.file has been sorted by p-values
# This algorithm does not allow the LD (not GWAS) snps presenting in mulitiple locus clusters
# Run scripts at /GWD/bioinfo/projects/statgen3/PGx6691_STOPGAP/STOPGAP2.3/scripts/locus.cluster
locus.cluster(ld.file="./LocusClustering/ld")

#Fine mapping of casual variants using PICS 
fine.map.pics(cluster.file = "./LocusClustering/STOPGAP_SNP_Cluster_",
gwas.file="stopgap_4sources_dbSNP141_msh.RData")

# For each unique GWAS SNP, rank its LD SNPs' genes based on their gene scores and r2 if gene scores tie
gene.rank(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_chr")

# Add cluser number data to the best LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_chr",
                        cluster.file = "./LocusClustering/STOPGAP_SNP_Cluster_")

#combine stopgap.data.bestld datasets
Combine.bestld.data(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster_chr")

# Generate the bestld data per gene-MeSH.
# Rank by the p-value bins and then by gene scores for  bestld data per gene-MeSH.
# P-value bins (-log10(P)): > 12 as 1; 8-12 as 2 and <8 as 3.
mesh.gene(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster.RData")

# merge.orphanet: Merge stopgap.gene.mesh with latest Orphanet data processed on 11/13/2104
# Create merge between Orphanet, and stopgap.gene.mesh 
# Change "Link" to "pubmedid" and add the OrphID as "pubmedid" when the Source=="Orphanet". 
# Set Rank=1, pval2 as 0 and GeneScore=max(stopgap.gene.mesh$gene.score)
# Remove MSH.TOP and OrphID columms and change the names to being consistent with stopgap.gene.mesh
merge.orphanet(gm.file = "./Mesh_gene/STOPGAP_r2_0.7_rmMHC_bestLD_gene_mesh.RData",
                       orphanet.file = "./Data/orphanet.RData",
					   gencode.file="./Data/gencode_v19_clean.RData")

revice.bestld(gm.file = "./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_MHC_genescore_generank_bestLD_cluster.RData",
                          gencode.file="./Data/gencode_v19_clean.RData" )


					   
