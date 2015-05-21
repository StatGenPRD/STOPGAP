#!/bin/bash
#Andrew Slater

CHROM=$1
BED=$2
OUTDIR=$3
LABEL_PY=Label_LD_results.py

#Run vcftools to calculate LD
zcat ALL.chr$CHROM.recode.vcf.gz | vcftools --vcf - --geno-r2-positions $BED --min-r2 0.5 --ld-window-bp 500000 --temp $OUTDIR --out $OUTDIR/chr$CHROM.LD

#Run python script to label LD
$LABEL_PY -d $OUTDIR -c $1
