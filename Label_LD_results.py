#!python2.7
from optparse import OptionParser
import logging
import re
import os
import pwd
import subprocess
import sys
import gzip

#Andrew Slater

#Directory containing the per-chromosome files with 1000 Genomes variant definitions (POS, ID, REF and ALT cut from VCF so order maintained)
#named chrN.ID.txt.gz
VARDEF = '.'

#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
#assumes files named chrN.LD.list.geno.ld and chrN.bed
parser.add_option('-d', '--ld-dir', help = 'Path to directory containing the LD results, will also be used to write labeled results',
                  action = 'store', type = 'string', dest = 'indir', default = '')
parser.add_option('-c', '--chromosome', help = 'Chromosome',
                  action = 'store', type = 'string', dest = 'chrom', default = '')

(options, args) = parser.parse_args()

#Get absolute paths of arguments
indir = os.path.realpath(options.indir)

#Check directory exist
if not os.path.isdir(indir) :
        raise Exception('Directory [ ' + indir + ' ] does not exist')

#Check LD results exist
if not os.path.isfile(indir + '/chr' + options.chrom + '.LD.list.geno.ld') :
        raise Exception('LD file [' + indir + '/chr' + options.chrom + '.LD.list.geno.ld] does not exist')


with gzip.open(VARDEF + '/chr' + options.chrom + '.ID.txt.gz','rb') as TG, open(indir + '/chr' + options.chrom + '.LD.list.geno.ld','r') as LD :
        TGVCF = {}
        #Single header CHR1 POS1 CHR2 POS2 N_INDV R^2
        LDi = LD.readline()
        #read through LD results and identify unique POS values that need labeling
        for LDi in LD :
                data = LDi.strip().split('\t')
                TGVCF[data[1]] = []
                TGVCF[data[3]] = []
        #Single header POS ID REF ALT
        vari = TG.readline()
        #read through ID file and store ID REF and ALT for each POS in LD results
        for c,vari in enumerate(TG) :
                coord = vari.strip().split('\t')
                if coord[0] in TGVCF.keys() :
                        TGVCF[coord[0]].append(coord)
 
with open(indir + '/chr' + options.chrom + '.LD.list.geno.ld','r') as LD , open(indir + '/chr' + options.chrom + '.labeled.ld','w') as out:
        #Initialize out header
        header = ['CHR1','POS1','CHR2','POS2','N_INDV','R^2','ID1','REF1','ALT1','ID2','REF2','ALT2']
        out.write('\t'.join(header) + '\n')
        #Single header CHR1 POS1 CHR2 POS2 N_INDV R^2
        LDi = LD.readline()
        #initialize tracker of prior result POS values
        P1 = 0
        P2 = 0
        #initialize indexes for determining replicant POS indicating variants with same POS (i.e. SNV and indel)
        I1 = 1
        I2 = 1
        #read through LD results and determine index based on previous result
        for LDi in LD :
                data = LDi.strip().split('\t')
                if int(data[1]) == P1 :
                        if int(data[3]) == P2  :
                                I2 += 1
                        else :
                                I2 = 1
                        if int(data[3]) < P2 :
                                I1 += 1
                else :
                        I1 = 1
                        I2 = 1
                #Store POS values for determining index of next result
                P1 = int(data[1])
                P2 = int(data[3])
                #Using index, lookup ID(s), REF and ALT for each POS1 and POS2
                ID1 = TGVCF[data[1]][I1 - 1][1].split(';')
                REF1 = TGVCF[data[1]][I1 - 1][2]
                ALT1 = TGVCF[data[1]][I1 - 1][3]
                ID2 = TGVCF[data[3]][I2 - 1][1].split(';')
                REF2 = TGVCF[data[3]][I2 - 1][2]
                ALT2 = TGVCF[data[3]][I2 - 1][3]
                #Restrict to results within 500 kb
                if abs(int(data[1]) - int(data[3])) <= 500000 :
                        #for each ID for POS1 and POS2, write input LD result with ID, REF and ALT
                        for i in ID1 :
                                for j in ID2 :
                                        datap = data + [i,REF1,ALT1,j,REF2,ALT2]
                                        out.write('\t'.join(datap) + '\n')
										