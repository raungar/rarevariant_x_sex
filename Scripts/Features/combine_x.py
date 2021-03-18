#!/bin/python

import argparse, re, gzip,glob, numpy
#from collections import defaultdict


parser=argparse.ArgumentParser()
parser.add_argument("--indir", type=str, help="the bed file", required=True)
parser.add_argument("--outfile_all", type=str, help="outfile for this sex containing bed info and zscores", required=True)
parser.add_argument("--outfile_inboth", type=str, help="outfile for this sex containing bed info and zscores", required=True)
parser.add_argument("--sexfile", type=str, help="the sex of the individual", required=True)
args = parser.parse_args()

outf_sex_all=gzip.open(args.outfile_all,"wb")
outf_sex_inboth=gzip.open(args.outfile_inboth,"wb")



#get sex per individual
sex_convert_key={}
sex_convert_key["1"]="male"
sex_convert_key["2"]="female"
sex_key={}
#get sex of individual
with open(args.sexfile,"r") as sex_f_read:
        sex_f_read.readline()
        for sex_line in sex_f_read.readlines():
                sex_line_split=sex_line.split("\t")
                sex_key[sex_line_split[0]]=sex_convert_key[sex_line_split[2]]

#compile at the beginning for speed
af_both=re.compile("AF_nfe=")

#key is position, value(m,f) -- if set len ==2 then seen in both sexes
seen_dic={}


####put positions into found dic
for this_f in glob.glob(args.indir+"/*bed"):
	print(this_f)
	with open(this_f, 'r') as f:
		f_split=this_f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=next(filter(lambda x: re.search('GTEX',x),f_split_again))
		sex=sex_key[ind]
		for line in f:
			line_split=line.split("\t")
			chr=line_split[0] #chr
			pos=line_split[1] #pos
			if pos in seen_dic:
				seen_dic[pos].add(sex)
			else:
				seen_dic[pos]={sex}
print(seen_dic)

for this_f in glob.glob(args.indir+"/*bed"):
	print(this_f)
	with open(this_f, 'r') as f:
		f_split=this_f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=next(filter(lambda x: re.search('GTEX',x),f_split_again))
		sex=sex_key[ind]

		print(ind)
		print(sex)
		for line in f:
			line_split=line.split("\t")
			chr=line_split[0] #chr
			pos=line_split[1] #pos
			gtex_maf=line_split[3] #GTEx MAF
			ensg=line_split[42]
			genetype=line_split[46]
			gnomad_split=(line_split[34]).split(';')
			if gnomad_split[0] == "NO_MATCH":
				print("NO MATCH: ",ref,",",alt,",",gtex_maf)
				#this prepares for cases where the reference allelse is actually teh minor allele
				if ref == alt:
					gnomad_maf_both=float(gtex_maf)
				#if it's been seen in gtex more than once, it shouldnt have MAF of zero					elif float(gtex_maf) > 0.001:
					gnomad_maf_both=float(gtex_maf)
				else:
					gnomad_maf_both=float(0)		
				print("both: ",gnomad_maf_both," , m: ", gnomad_maf_m, " , f: ", gnomad_maf_f)
			else:
				gnomad_maf_both=float((([col for col in gnomad_split if af_both.match(col)])[0].split("="))[1])
			#print("\t".join([chr,pos,ensg, genetype,ind,sex,gtex_maf,str(gnomad_maf_both)]))
			myline=[chr,pos,pos,str(gnomad_maf_both),ind,"SNPs",ensg, genetype,sex,gtex_maf]
			#must be seen in both
			if (len(seen_dic[pos])<2):
				continue
				outf_sex_all.write(("\t".join(myline)+"\n").encode('utf-8'))
			else:
				outf_sex_inboth.write(("\t".join(myline)+"\n").encode('utf-8'))
outf_sex_all.close()
outf_sex_inboth.close()
