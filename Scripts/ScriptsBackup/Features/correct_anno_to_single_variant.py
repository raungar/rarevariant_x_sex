##########purpose: take only smallest maf and get gene-level ino
##########input: ind maf files, params of X percent difference of MAF F/M
##########output: 
############(1) file of MAF difference between M/F of more than X percent
############collapsed files for (2) males (3) females and (4) both that is as follows:
############ when referring to mafs, unless explicitly otherwise this refers to the min MAF
############gene/ count of MAF<0.01 for ind /  MAF GTEX/ MAF GNOMAD nonFIN BOTH / 
############MAF GNOMAD nonFIN F/ MAF GNOMAD nonFIN M / minChr / Min Pos /
############ minRef / MinALT / minRefALT / min percent differ / geneTYPE
import argparse
import glob
import re
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--infile', type=str, help='file  to read from')
parser.add_argument('--outfile', type=str, help='file to write to')
parser.add_argument('--chr', type=str, help='chr number, write x if x chromosome :) ')

args = parser.parse_args()
infile=args.infile
outfile=args.outfile
this_chr=args.chr
outwrite=gzip.open(outfile,"wb")

#compile at the beginning for speed
af_gnomad=re.compile("AF_nfe=")

print(infile)
#open file and read line by line
#get individual ID and sex of this individual
#key is var, value is max maf
dic_maf=dict()
dic_all=dict() #key is var, value is line

#so can use with both x and aut
if str(this_chr).lower() == "x":
	gnomad_anno_pos=34
else:
	gnomad_anno_pos=42

#actually open file, read line by line
with open(infile,"r") as f_read:
	for line in f_read.readlines():
		line_split=line.split("\t")
		#get interesting columns
		chr=line_split[0] #chr
		pos=line_split[1] #pos
		gtex_maf=line_split[3] #GTEx MAF
		ref=line_split[5] #reference allele
		alt=line_split[6] #alternate allele
		#if(line_split[42] == "NA"):
		#	next
		gnomad_anno=line_split[gnomad_anno_pos] #42] # GNOMAD annotation to further split for AUT
		#gnomad_anno=line_split[] #34] # GNOMAD annotation to further split for X chr
		gnomad_split=gnomad_anno.split(';')
		if gnomad_split[0] == "NO_MATCH":
			#this prepares for cases where the reference allelse is actually teh minor allele
			if ref == alt:
				gnomad_maf=float(gtex_maf)
			#if it's been seen in gtex more than once, it shouldnt have MAF of zero
			elif float(gtex_maf) >= 0.01:
				gnomad_maf=float(gtex_maf)
			else:
				gnomad_maf=float(0)		
		else:
			found_af=[col for col in gnomad_split if af_gnomad.match(col)]
			try:
				gnomad_maf=float(((found_af[0]).split("="))[1])
			except:
				gnomad_maf=float(gtex_maf)
			##gnomad_maf=float(([0].split("="))[1])
			#gnomad_maf=float((([col for col in gnomad_split if af_gnomad.match(col)])[0].split("="))[1])

		#if there is already a RV recorded for this gene
		if pos in dic_maf:
			#get current max MAF
			dic_current_max_maf=float(dic_maf[pos])
			#if this is more common, store this instead
			if gnomad_maf > dic_current_max_maf:
				dic_maf[pos]=dic_current_max_maf
				dic_all[pos]=line_split
		#if pos not seen yet
		else:
			dic_maf[pos]=float(gnomad_maf)
			dic_all[pos]=line_split

##writes the dictionary to a file
for pos_key in dic_all:	
	outwrite.write(('\t'.join(map(str,dic_all[pos_key]))).encode())
outwrite.close()
