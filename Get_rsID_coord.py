#!python2.7
#Andrew Slater
from optparse import OptionParser
import logging
import re
import os
import pwd
import subprocess
import sys



#extract arguments
parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-r', '--rs-list', help = 'Path to text file listing the rs IDs with optional #-prefixed comment lines',
                  action = 'store', type = 'string', dest = 'rslist', default = '')
parser.add_option('-o', '--output-directory', help = 'Path to an existing output directory where the results will be written, any existing files with same name will be overwritten',
                  action = 'store', type = 'string', dest = 'outdir', default = '')

(options, args) = parser.parse_args()

#Get absolute paths of arguments
outdir = os.path.abspath(options.outdir)
rslistf = os.path.abspath(options.rslist)

#Create output directory
if not os.path.isdir(outdir):
	try:
		os.mkdir(outdir)
	except OSError:
		raise Exception('Directory [ ' + outdir + ' ] does not exist and cannot be created')


#create log
logging.basicConfig(filename=outdir + '/Get_rsID_coord.log',filemode='w', format='%(levelname)s\t%(asctime)s\t%(message)s', level=logging.DEBUG)

#start log with user and arguments - will validate file contents as process
logging.info('user [%s;%s] looking up GRCh37 coordinates from [%s] and reporting in [%s]', pwd.getpwuid(os.getuid())[0], pwd.getpwuid(os.getuid())[4], rslistf, outdir)

#take list of rs IDs and standardize to lowercase rs followed by digits without any leading zeros
rscoord = {}
rsIDmap = {}
count = 0
notrscount = 0
with open(rslistf, "r") as rslist, open(outdir + '/notrsIDs.txt', "w") as notrs :
	notrs.write('value_from_input_rs_list\tline_num_from_input_rs_list\n');
	for line_num, rsid in enumerate(rslist) :
		rsid = rsid.strip()
		if rsid == '' or rsid[0] == '#' :
			continue
		elif re.search(r"^rs\d+$", rsid.lower()) == None :
                        #Produces too many lines in log file which is redundant with the output file, instead reporting single count below
			#logging.warning('Line number [%d] not a valid rs ID [%s] writing to [%s]', line_num, rsid, outdir + '/notrsIDs.txt')
                        notrscount += 1
			notrs.write(rsid + '\t' + str(line_num) + '\n')
		else :
			id = re.search(r"^rs(\d+)$", rsid.lower())
			#remove leading zeros
			intid = int(id.group(1))
			strid = str(intid)
			#initialize coordinate to empty
			rscoord[strid] = ''
			#store all unique input IDs by standardized ID
			#If standardized ID already exists but this is first occurance of input ID, then append
			if 'rs' + strid in rsIDmap and rsid not in rsIDmap['rs' + strid] :
				count += 1
				rsIDmap[strid].append(rsid)
			#If this is first occurance of standardized ID, record input ID
			elif 'rs' + strid not in rsIDmap :
				count += 1
				rsIDmap[strid] = [rsid]
#
logging.warning('Found [%d] invalid input values where no rs ID could be extracted, written to [%s]', notrscount, outdir + '/notrsIDs.txt')
	
logging.info('Found [%d] unique valid input values corresponding to [%d] unique rs IDs in [%s]', count, len(rscoord), rslistf)

#HGVS chromosomes for GRCh37:
chroms = {"NC_000001.10":"1",
        "NC_000002.11":"2",
        "NC_000003.11":"3",
        "NC_000004.11":"4",
        "NC_000005.9":"5",
        "NC_000006.11":"6",
        "NC_000007.13":"7",
        "NC_000008.10":"8",
        "NC_000009.11":"9",
        "NC_000010.10":"10",
        "NC_000011.9":"11",
        "NC_000012.11":"12",
        "NC_000013.10":"13",
        "NC_000014.8":"14",
        "NC_000015.9":"15",
        "NC_000016.9":"16",
        "NC_000017.10":"17",
        "NC_000018.9":"18",
        "NC_000019.9":"19",
        "NC_000020.10":"20",
        "NC_000021.8":"21",
        "NC_000022.10":"22",
        "NC_000023.10":"X",
        "NC_000024.9":"Y"}

#looping through sets of 100 for Entrez efetch, but if a single id is invalid, the entire fetch would fail
#so need to check ids first
i = 0
while i < len(rscoord) :
        j = min(i + 100, len(rscoord))
        cmd = "efetch -db snp -id {} -mode asn.1 | grep 'rsId' | sed -e 's/^\s*rsId\s*//g' | sed -e 's/,//' > {}".format(','.join(rscoord.keys()[i:j]), outdir + '/fetch.out')
        #print cmd
        os.system(cmd)
        fetch = open(outdir + '/fetch.out', 'r')
        idfound = []
        for line in fetch :
                line = line.strip()
                if line in rscoord.keys()[i:j] :
                        idfound.append(line)
        fetch.close()
        cmd = "efetch -db snp -id {} -mode XML -format docsum | egrep '<DocumentSummary><Id>|<SNP_ID>|<DOCSUM>' > {}".format(','.join(idfound), outdir + '/fetch.out')
        #print cmd
        os.system(cmd)
        fetch = open(outdir + '/fetch.out', 'r')
        for line in fetch :
                if line[0:21] == '<DocumentSummary><Id>' :
                        line = line.split('<')
                        rsid = line[2][3:]
                        #print rsid
                elif line[0:9] == '\t<SNP_ID>' :
                        line = line.split('<')
                        currid = line[1][7:]
                elif line[0:9] == '\t<DOCSUM>' :
                        #print line
                        line = line.split('|')
                        for m in line :
                                parse1 = re.search(r"HGVS=(.+)$", m)
                                if parse1 != None :
                                        line2 = parse1.group(1).split(',')
                                        #print line2
                                        found = '0'
                                        for k in line2 :
                                                #print "k:" + k
                                                #print "found:" + "|".join(found)
                                                k = k.split(':')
                                                if k[0] in chroms :
                                                  if len(k) > 1 and k[1][0:2] == 'g.' :
                                                        chrom = chroms[k[0]]
                                                        parse = re.search(r"^g\.(\d+)(.+)>(.+)$", k[1])
                                                        if parse == None :
                                                                vtype = 'not_snv'
                                                                coord = ''
                                                                ref = ''
                                                                alt = ''
                                                        else :
                                                                coord = parse.group(1)
                        #                                             print coord
                                                                ref = parse.group(2)
                                                                alt = parse.group(3)
                        #                                           print ref
                        #                                            print alt
                                                                if len(ref) > 1 or len(alt) > 1 :
                                                                        vtype = 'not_snv'
                                                                else :
                                                                        vtype = 'snv'
                                                        #print chrom + ':' + coord + vtype + ref + '|' + alt
                                                        #print found
                                                       #if no coordinates found yet, initialize to this coord
                                                        if found == '0' :
                                                                #print "in found0"
                                                                found = [chrom + ':' + coord, vtype, ref, alt, 'rs' + currid]
                                                        #if already set to blank due to observing more than one coordinate, skip to next
                                                        elif found[0] == '' :
                                                                #print "in found1"
                                                                continue
                                                        #if in PAR and Y coord found but now have X, replace with X
                                                        elif found[0][0] == 'Y' and coord != '' and ((int(coord) >= 60001 and int(coord) <= 2699520) or (int(coord) >= 154931044 and int(coord) <= 155260560 )) and chrom == 'X' :
                                                                #print 'Xhere'
                                                                found = [chrom + ':' + coord, vtype, ref, alt, 'rs' + currid]
                                                        #if in PAR and X coord found but now have Y, skip to next
                                                        elif found[0][0] == 'X' and coord != '' and ((int(coord) >= 10001 and int(coord) <= 2649520) or (int(coord) >= 59034050 and int(coord) <= 59363566 )) and chrom == 'Y' :
                                                                #print 'here'
                                                                continue
                                                        #if more than one coordinate, set blanks
                                                        else :
                                                                #print "found5"
                                                                found = ['', vtype, '', '', 'rs' + currid]
                                                  else :
                                                          #print "found6"
                                                          found = ['','','','','rs' + currid]
                                        if found != '0' :                                
                                                rscoord[rsid] = found
                                                #print 'final: ' + str(rsid) + ' ' + str(found)
                          
        fetch.close()
        i = j
#print(rscoord)
os.remove(outdir + '/fetch.out')

#determine what not found and write output file
bed = {}
with open(outdir + '/rsID_Coordinates.txt',"w") as outlist, open(outdir + '/rsIDsNotFound.txt',"w") as noutlist, open(outdir + '/rsIDNotSNV.txt', "w") as ioutlist, open(outdir + '/rsIDambig.txt', "w") as aoutlist :
	ncount = 0
	fcount = 0
	icount = 0
	acount = 0
	outlist.write('#rsID_searched\tGRCh37_Coordinate\ttype\tREF\tALT\tcurrent_rsID\tinputrsIDs\n')
	noutlist.write('#rsID_searched\tinputrsIDs\n')
	ioutlist.write('#rsID_searched\tGRCh37_Chromosome\ttype\tcurrent_rsID\tinputrsIDs\n')
	aoutlist.write('#rsID_searched\ttype\tcurrent_rsID\tinputrsIDs\n')
	for rs in rscoord.keys() :
		if rscoord[rs] == '' :
			ncount += 1
			noutlist.write('rs' + rs + '\t' + ";".join(rsIDmap[rs]) + '\n')
		elif rscoord[rs][1] == 'not_snv' :
                        icount += 1
                        oline = ['rs' + rs, rscoord[rs][0], rscoord[rs][1], rscoord[rs][4], ";".join(rsIDmap[rs])]
                        ioutlist.write('\t'.join(oline) + '\n')
			
                        #dbSNP VCFs report PAR on chr Y, converting to chrX as 1000 Genomes data
			#PAR = re.search(r"^Y:(\d+)",rscoord[rs][0])
			#if PAR != None :
			#	Ycoord = int(PAR.group(1))
			#	if Ycoord >= 10001 and Ycoord <= 2649520 :
			#		rscoord[rs][0] = 'X:' + str(Ycoord + 50000)
			#	elif Ycoord >= 59034050 and Ycoord <= 59363566 :
			#		rscoord[rs][0] = 'X:' + str(Ycoord + 95896994)
			
		elif rscoord[rs][0] == '' :
                        acount += 1
                        oline = ['rs' + rs, rscoord[rs][1], rscoord[rs][4], ';'.join(rsIDmap[rs])]
                        aoutlist.write('\t'.join(oline) + '\n')
                else :
                        fcount += 1
                        outlist.write('rs' + rs + '\t' + "\t".join(rscoord[rs]) + '\t' + ";".join(rsIDmap[rs]) + '\n')
                        #add coordinate to bed dictionary
                        chrpos = rscoord[rs][0].split(':')
                        #If first coordinate for this chromosome then initialize list of positions
                        if chrpos[0] not in bed :
                                bed[chrpos[0]] = [int(chrpos[1])]
                        #Append this position to the existing list
                        else :
                              bed[chrpos[0]].append(int(chrpos[1]))

#Create bed files for each chromosome
#Note, for vcftools v0.1.12b --geno-r2-positions, start and end must both = POS instead of traditional bed format where start is POS - 1
for chrom in bed.keys() :
        with open(outdir + '/chr' + chrom + '.bed',"w") as bedout :
                line = ('#chrom','start','end')
                bedout.write('\t'.join(line) + '\n')
                bed[chrom].sort()
                for coord in bed[chrom] :
                        line = (chrom,str(coord),str(coord))
                        bedout.write('\t'.join(line) + '\n')

if fcount > 0 :
	logging.info('Found coordinates for a total of [%d] SNV rs IDs, written to [%s]', fcount, outdir + '/rsID_Coordinates.txt')
else :
        os.remove(outdir + '/rsID_Coordinates.txt')
if icount > 0 :
        logging.warning('Found [%d] non-SNV rs IDs where the VCF coordinates and alleles could not be identified from the HGVS notation, written to [%s]', icount, outdir + '/rsIDNotSNV.txt')
        logging.warning('You may be able to look these up manually to determine correct coordinate or alleles and then add to [%s]', outdir + '/rsID_Coordinates.txt')
        logging.warning('Note, if the GRCh37_Chromosome is blank it is because this non-SNV variant is also ambiguously mapped or has more than one alternative allele')        
else :
        os.remove(outdir + '/rsIDNotSNV.txt')
if acount > 0 :
        logging.warning('Found [%d] rs IDs mapping to more than one coordinate or more than one alternative allele, written to [%s]', acount, outdir + '/rsIDambig.txt')
        logging.warning('You may be able to look these up manually to determine correct coordinate or allele and then add to [%s]', outdir + '/rsID_Coordinates.txt')
        logging.warning('Note, the chr Y coordinate for variants in the pseudo-autosomal region of chr X are ignored and thus considered uniquely mapped')
        
else :
	os.remove(outdir + '/rsIDambig.txt')
if ncount > 0 :
	logging.warning('Did not find [%d] rs IDs in dbSNP via EntrezDirect, writing to [%s]', ncount, outdir + '/rsIDsNotFound.txt')
	logging.warning('You may be able to look these up manually (usually they have been merged into another rs ID) and then add to [%s]', outdir + '/rsID_Coordinates.txt')
else :
	os.remove(outdir + '/rsIDsNotFound.txt')
if len(bed) > 0 :
        logging.info('Created bed files for [%d] chromosomes', len(bed))
        
