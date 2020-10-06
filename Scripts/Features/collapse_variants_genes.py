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

#parser = argparse.ArgumentParser(description='argparser')
#parser.add_argument('dir_read', type=str, help='directory to read file from')
#parser.add_argument('out_mafdiff', type=str, help='Output file for MAF Mvs.F diff')
#parser.add_argument('cutoff_mafdiff', type=float, help='cutoff [0-1] for % difference')
#parser.add_argument('out_m', type=str, help='outfile collaped m')
#parser.add_argument('out_f', type=str, help='outfile collaped f')
#parser.add_argument('out_b', type=str, help='outfile collaped both')
#args = parser.parse_args()

dir_read="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/bySiteAnnoX/GenesAnno"

#compile at the beginning for speed
af_both=re.compile("AF_nfe=")
af_m=re.compile("AF_nfe_male=")
af_f=re.compile("AF_nfe_female=")

### loop through SNPS/INDELS/SV
for vartype in ["SNP","indel","SV"]:
	## loop through all inds in this vartype
	for f in glob.glob(dir_read+"/*"+vartype+"*"):
		#open file and read line by line
		f_split=f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=f_split_again[1]
		with open(f,"r") as f_read:
			for line in f_read.readlines():
				line_split=line.split("\t")
				chr=line_split[0] #chr
				pos=line_split[1] #pos
				gtex_maf=line_split[3] #GTEx MAF
				ref=line_split[5] #reference allele
				alt=line_split[6] #alternate allele
				ensg=(line_split[42]).split("\\"") # ENSG ID
				genetype=line_split[46] #lincRNA, proteincoding, etc
				varswitch=line_split[53].strip() ## A->C as AC
				gnomad_anno=line_split[34] # GNOMAD annotation to further split
				gnomad_split=gnomad_anno.split(';')
				gnomad_maf_both=float((([col for col in gnomad_split if af_both.match(col)])[0].split("="))[1])
				gnomad_maf_m=float((([col for col in gnomad_split if af_m.match(col)])[0].split("="))[1])
				gnomad_maf_f=float((([col for col in gnomad_split if af_f.match(col)])[0].split("="))[1])
	
				gnomad_maf_diff=(gnomad_maf_m-gnomad_maf_f)/((gnomad_maf_m+gnomad_maf_f)/2)

				#if gnomad_maf_diff > .1:
					#write to gnomad diff file this information
				print("\t".join([chr,pos,ensg, vartype,ind,ref,alt,varswitch,
						gtex_maf, str(gnomad_maf_both),str(gnomad_maf_m),str(gnomad_maf_f),
						genetype,str(gnomad_maf_diff)]))

				




