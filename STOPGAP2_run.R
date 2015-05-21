#--------Run STOPGAP2 R functions to generate STOPGAP2 datasets -----------------#
# Author: Judong Shen

setwd("STOPGAP.working.directory")
library(gplots)
library(plyr)
library(reshape)
#library(cba)
#library(ggplot2)

# (1) merge the grasp, nhgri and gwasdb datasets
p.threshold <- 1E-04
grasp <- grasp.import(file = "GraspFullDataset2.zip", p.thres=p.threshold)
nhgri <- NHGRI.import(file = "gwascatalog.txt", p.thres=p.threshold)
gwasdb <- gwasdb.import(file = "gwasdb_20140812_snp_trait.gz", p.thres=p.threshold)
data <- data.merge(nhgri, grasp, gwasdb)

gwas.gap(new.gwas="stopgap_3sources.RData", old.gwas="../../STOPGAP2_pipeline/Data/stopgap_3sources.RData")

## Manually generate the disease.msh.txt file which includes all the previous msh/category information and the fuzzy matched new information
## Manually fill the gap and then generate the final table called disease.msh.final.txt (for STOPGAP2 data)

# Identify the unique GWAS SNPs 
gwas.snps(gwas.file="stopgap_3sources.RData")

# Coordinate lookup and LD calculation based on the 1KG data for the new SNPs identified in this version
# Run LD under: /GWD/bioinfo/databases/popgen/hla/LD instead of GAF. 
run.ld()

# Update rsIDs in the gwas data (data) to dbSNP141 version
rsID.update(gwas.file="stopgap_3sources.RData", 
                        coor.file="./STOPGAP2_LDResults/rsID_Coordinates.txt")

#Process ChIA-PET files and merge them together
ChIA.PET.import(ChIA.PET.path="./ChIA-PET", 
                gencode.file="gencode_v20.RData",
                win.size=1000)

# Find all variants in LD with GWAS variants at r2 > 0.5 - GWAS, previous version of LD data, the path to the new additional LD data
ld.snps(gwas.data="stopgap_3sources_dbSNP141.RData",
                    ld.path="./STOPGAP2_LDResults")
  
# updata.ld.data: updata ld data:
# Remove SNP.ld with two positions in both ld.snps and ld.snps.r2
# Use chr:pos to replace the missing ld SNP ids in both ld.snps and ld.snps.r2
updata.ld.data(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.RData", 
                ld.r2.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.RData")


# Update the ld.snps.r2 and ld_snps data by including the GWAS SNPs data which don't have any LD information there
update.ld.snps.r2(gwas.file="stopgap_3sources_dbSNP141.RData",
                              ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.RData",
                              rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt")

# Update the ld.snps.r2 and ld_snps data by removing the LD SNPs more than 500kb away from the GWAS SNPs,
# If any GWAS SNPs get deleted, put them back by including r2=1 there
update1.ld.snps.r2(gwas.file="stopgap_3sources_dbSNP141.RData",
                               ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.RData",
                               dis = 500*1000)

# Update the ld.snps.r2 data by adding the AF information from 1KG phase I data ...
r2.frq(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.RData",
                   frq.file="./1KG_AF/gws.frq.RData")

# ld.rdb.dhscor: add the Cat.rdb to the ld snps &
# For each LD SNP, (1) identify each dhscor where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.rdb.dhscor(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
                          rdb.file="./ENCODE/regulomedb14.RData",
                          dhscor.file="./ENCODE/dhscor.RData",
                          mpos.rng.p = c("Chr.pdhs", "Start.pdhs", "Stop.pdhs"),
                          mpos.rng.d = c("Chr.ddhs", "Start.ddhs", "Stop.ddhs"))

# ----- Merge the CHIA-PET data
# ld.chiapet: For each LD SNP, (1) indentify each chiapet where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.chiapet(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
           chiapet.file="./CHIA_PET/CHIA.PET.RData",
           mpos.rng.p = c("pChr", "pStart", "pEnd"),
           mpos.rng.d = c("dChr", "dStart", "dEnd"))

# ld.eQTL: For each LD SNP, (1) indentify each eQTL where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, cor (if in promoter, set cor=1), and an indicator if it is in a distal or promoter region.
ld.eQTL(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
                    eQTL.file="./eQTL/eQTL_UChichago/eQTL_UChicago_STOPGAP.RData",
                    eQTL.GTEx.file="./eQTL/GTEx/GTEx.data.RData",
                    eQTL.GTEx.sz="./eQTL/GTEx/GETx_Tissue_SampleSize.txt")


# ld.fantom5: For each LD SNP, (1) indentify each fantom5 where SNP falls within promoter or distal genomic region; 
#                             (2) add gene, tissue and an indicator if it is in a distal region (note  no promoter region here).
ld.fantom5(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData", 
           fantom5.file="./FANTOM5/FANTOM5.RData",
           mpos.rng.d = c("chr", "startpos", "endpos"))

# ld.Posgene: For each LD SNP, (1) identify genes falling in the 5Kb windows (upstream 5kb, downstream 5kb) 
ld.Posgene(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
           Posgene.file="./gencode_v20_RefSeq37_1.RData",
           win.size = 5000)

# ld.vep: Merge gwas ld SNPs and VEP information. 
ld.vep(ld.file="./STOPGAP2_LDResults/gwas_LD/ld.snps.updated.plusGWASsnps.500kb.RData",
       vep.path="./1KG/Annotations")

# ld.vep.simplify: Simplify the full version of merged gwas ld SNPs and VEP data. 
ld.vep.simplify(vep.file="./VarGeneMapping/gwas.ld.vep.full.RData", 
                severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")

#merge.ld.simp: Merge based on SNP.ld and gene from each data
# % of unique (ld) GWAS SNPs with at least one genes
merge.ld.simp(dhscor.file ="./VarGeneMapping/gwas.ld.rdb.dhscor.RData",
              chiapet.file = "./VarGeneMapping/gwas.ld.chiapet.RData",
              fantom5.file = "./VarGeneMapping/gwas.ld.fantom5.RData",
              eQTL.file = "./VarGeneMapping/gwas.ld.eQTL.RData",
              posgene.file = "./VarGeneMapping/gwas.ld.PosGene.RData",
              vep.file.simplified = "./VarGeneMapping/gwas.ld.vep.simplify.RData")

# Subset of the var2gene data by using the union of RefSeq and Gencode gene names
# Also Update the var2gene mapping data by removing any rows with conflicting chr.ld with those from the chr information from gencode data
var2gene.filter(var2gene.file = "./VarGeneMapping/gwas.ld.var2gene.simplified.RData",
                gencode.file = "./gencode_v20.RData",
                RefSeq.file = "./RefSeq37_1.RData")

# gwas.mesh: add mesh terms to the GWAS data. 
# Run it on Windows only
gwas.mesh(gwas.file="./stopgap_3sources_dbSNP141.RData",
          disease2mesh.file = "./disease_mesh.txt")

# Merge GWAS data, LD r2 data and var2gene_vepSimp data together
gwas.ld.var2gene1(gwas.file="stopgap_3sources_dbSNP141_msh.RData",
                  ld.file = "./STOPGAP2_LDResults/gwas_LD/ld.snps.r2.updated.plusGWASsnps.500kb.AF.RData",
                  var2gene.file = "./VarGeneMapping/Var2Gene_vepSimp.RData")

rm.mhc.genes(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7.RData",
             lRmMHC=TRUE)

# Add gene scores to the final Mergedata_VEPsimplified_withGene_r2_0.7.RData data
gene.score(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC.RData",
           severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")

# For each unique GWAS SNP, rank its LD SNPs' genes based on their gene scores and r2 if gene scores tie
gene.rank(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore.RData")

# gene.score1(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5_rmMHC.RData",
#             severity.file = "./VarGeneMapping/VEP_Consequence_w_severity.txt")
# 
# gene.rank1(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.5_rmMHC_genescore.RData")

#ld.data: Get the LD data as a starting point for SNP locus clustering:
ld.data(data.file="./LocusClustering/gwas_r2.RData")
locus.cluster1(ld.file="./LocusClustering/ld1.RData")
locus.cluster2(ld.file="./LocusClustering/ld2.RData")
locus.cluster3(ld.file="./LocusClustering/ld3.RData")
locus.cluster4(ld.file="./LocusClustering/ld4.RData")
locus.cluster5(ld.file="./LocusClustering/ld5.RData")
cluster.merge(c.file1="./LocusClustering/ld1_LocusClusters.RData",
              c.file2="./LocusClustering/ld2_LocusClusters.RData",
              c.file3="./LocusClustering/ld3_LocusClusters.RData",
              c.file4="./LocusClustering/ld4_LocusClusters.RData",
              c.file5="./LocusClustering/ld5_LocusClusters.RData")

# Add cluser number data to the full LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank.RData",
            cluster.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData")

# Add cluser number data to the best LD data
# Also add the ranks based on the p-values of snp.gwas within each cluster in each pubmedid/disease combination
add.cluster.bestld(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD.RData",
            cluster.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData")

# Generate the bestld data per gene-MeSH.
# Rank by the p-value bins and then by gene scores for  bestld data per gene-MeSH.
# P-value bins (-log10(P)): > 12 as 1; 8-12 as 2 and <8 as 3.
mesh.gene(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData")

# Reorder all related files based on columns
# Reorder for the bestLD version data
reorder.data(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_bestLD_cluster.RData")

# Reorder for the full version data
reorder.data1(data.file="./GWAS_LD_var2gene/Mergedata_VEPsimplified_withGene_r2_0.7_rmMHC_genescore_generank_cluster.RData")

# Rename the R object names by using simpler, descriptive names.
rename.data()

# Update the snp.cluster data by only keeping three columns "snp.ld", "cluster", "start.gwassnp"
snp.cluster.update(data.file = "./LocusClustering/STOPGAP_SNP_Clusters.RData")

# Update the gwas data by adding the chr and pos information
gwas.update(gwas.file = "../STOPGAP2/gwas.RData",
            rsid.file = "./STOPGAP2_LDResults/rsID_Coordinates.txt")

# Update all related STOPGAP2 results
results.update(res.path="../STOPGAP2")

# merge.omim1: Merge stopgap.gene.mesh with latest OMIM/Orphanet data processed on 11/13/2104
# Create merge between OMIM, Orphanet, and stopgap.gene.mesh -> stopgap.gwas.omim
# Change "Link" to "pubmedid" and add the OrphID as "pubmedid" when the Source=="Orphanet". 
# Set Rank=1, pval2 as 0 and GeneScore=max(stopgap.gene.mesh$gene.score)
# Remove MSH.TOP and OrphID columms and change the names to being consistent with stopgap.gene.mesh
# Merge the two datasets by using all column names of omim.nomhc data. 
#  Basically it is adding all rows of omim.nomhc to the stopgap.gene.mesh data 
merge.omim1(gm.file = "../STOPGAP2/stopgap.gene.mesh.RData",
            omim.file = "./OMIM_Orphanet/omim.RData")


# Remove the pubmedid=="22377632" or msh=="transmission" in the stopgap.RData, stopgap.bestld.RData and stopgap.gene.mesh.RData
# #Removing PUBMEDID from GWASDB only: PUBMEDID 21150878 & PUBMEDID 22190364
res.update1(rm.pubmedid="22377632",
            rm.pubmedid1=c("21150878","22190364"))

# Further add the function genomics and annotation evidence column to the bestld or gene.mesh data 
load("../STOPGAP2/stopgap.bestld.RData")
stopgap.bestld <- add.evid(stopgap.bestld)
names(stopgap.bestld) <- c("gene","msh","msh.tree","msh.cat","disease",
                           "pvalue","gene.score","gene.rank.max",
                           "gene.rank.min","gwassnp.p.rank","snp.ld", "snp.gwas",
                           "r2","cluster","clust.init.gwassnp",
                           "chr.ld","pos.ld","snp.gwas.orig",
                           "pubmedid","init.samp","source",
                           "rep.samp","nhlbikey","loc.paper",
                           "datepub","gwas.ancestry","sampsize.dis.rep",
                           "sampsize.dis","sampsize.rep","snp.gwasdb",
                           "source.gwasdb","hpo.id","hpo.term",
                           "do.id","do.term","af.1kg",
                           "asn.af.1kg","amr.af.1kg","afr.af.1kg",
                           "eur.af.1kg","cat.rdb","dhs.cor",
                           "dhs.type","chiapet.cell","chiapet.type",
                           "fantom5.tissue","fantom5.type","eqtl.ref",
                           "eqtl.tissue","eqtl.scoretype","eqtl.score",
                           "vep.canonical","vep.conseq","vep.aa",
                           "vep.condel.pred","vep.condel.score","vep.severity",
                           "evidence","gene.best","gene.score.best",
                           "evidence.best")
save(stopgap.bestld, file = "../STOPGAP2/stopgap.bestld.RData")
write.table(stopgap.bestld,
            "../STOPGAP2/stopgap.bestld.txt",
            sep = "\t", row.names = FALSE, quote = F, na = "")    


load("../STOPGAP2/stopgap.gene.mesh.RData")
stopgap.gene.mesh <- add.evid(stopgap.gene.mesh)
names(stopgap.gene.mesh) <- c("gene","msh","msh.tree","msh.cat","disease",
                              "pvalue","gene.score","gene.rank.max",
                              "gene.rank.min","gwassnp.p.rank","snp.ld", "snp.gwas",
                              "r2","cluster","clust.init.gwassnp",
                              "chr.ld","pos.ld","snp.gwas.orig",
                              "pubmedid","init.samp","source",
                              "rep.samp","nhlbikey","loc.paper",
                              "datepub","gwas.ancestry","sampsize.dis.rep",
                              "sampsize.dis","sampsize.rep","snp.gwasdb",
                              "source.gwasdb","hpo.id","hpo.term",
                              "do.id","do.term","af.1kg",
                              "asn.af.1kg","amr.af.1kg","afr.af.1kg",
                              "eur.af.1kg","cat.rdb","dhs.cor",
                              "dhs.type","chiapet.cell","chiapet.type",
                              "fantom5.tissue","fantom5.type","eqtl.ref",
                              "eqtl.tissue","eqtl.scoretype","eqtl.score",
                              "vep.canonical","vep.conseq","vep.aa",
                              "vep.condel.pred","vep.condel.score","vep.severity",
                              "asso.count","num.cluster","num.gwas.snp",
                              "num.ld.snp","num.gene.cluster","max.gene.score",
                              "evidence","gene.best","gene.score.best",
                              "evidence.best")
save(stopgap.gene.mesh, file = "../STOPGAP2/stopgap.gene.mesh.RData")
write.table(stopgap.gene.mesh,
            "../STOPGAP2/stopgap.gene.mesh.txt",
            sep = "\t", row.names = FALSE, quote = F, na = "")    

# row and col numbers of STOPGAP2 main datasets.
row.col()
