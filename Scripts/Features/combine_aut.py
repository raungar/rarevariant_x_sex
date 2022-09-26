#!/bin/python

import argparse, re, gzip,glob, numpy, os
#from collections import defaultdict


parser=argparse.ArgumentParser()
parser.add_argument("--indir", type=str, help="the bed file", required=True)
parser.add_argument("--genewindow", type=str, help="distance from gene start site to include", required=True)
parser.add_argument("--outfile_all", type=str, help="outfile for this sex containing bed info and zscores", required=True)
parser.add_argument("--outfile_inboth", type=str, help="outfile for this sex containing bed info and zscores", required=True)
parser.add_argument("--sexfile", type=str, help="the sex of the individual", required=True)
parser.add_argument("--gtf", type=str, help="gtf file for exon adding. assuming sorted for speed.", required=True)
parser.add_argument("--filename_match", type=str, help="keyword to find right file", required=True)
args = parser.parse_args()

outf_sex_all=gzip.open(args.outfile_all,"wb")
outf_sex_inboth=gzip.open(args.outfile_inboth,"wb")

genewindow=int(args.genewindow)

#open gtf fiel for reading, reading line by line (sorting is assumed for speed)
gtf=open(args.gtf,"r")
gtf_line=(gtf.readline()).split("\t")
gtf_start=int(gtf_line[3])

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
variant_type=re.compile("variant_type=")

#key is position, value(m,f) -- if set len ==2 then seen in both sexes
seen_dic={}


####put positions into found dic
for this_f in glob.glob(args.indir+"/*"+args.filename_match):
	print(this_f)
	with gzip.open(this_f, 'r') as f:
		f_split=this_f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=next(filter(lambda x: re.search('GTEX',x),f_split_again))
		sex=sex_key[ind]
		for line in f:
			line_split=line.decode('utf-8').split("\t")
			chr=line_split[0] #chr
			pos=line_split[1] #pos
			if pos in seen_dic:
				seen_dic[pos].add(sex)
			else:
				seen_dic[pos]={sex}

for this_f in glob.glob(args.indir+"/*"+args.filename_match):
	print(this_f)
	if (os.path.isdir(this_f)):
		continue
	with gzip.open(this_f, 'r') as f:
		f_split=this_f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=next(filter(lambda x: re.search('GTEX',x),f_split_again))
		sex=sex_key[ind]

		print(ind)
		print(sex)
		for line in f:
			line_split=line.decode('utf-8').split("\t")
			#print(line_split)
			chr=line_split[0] #chr
			pos=line_split[1] #pos
			gtex_maf=line_split[3] #GTEx MAF
			geno=line_split[4] #0/1/2 genotype
			ref=line_split[5]
			alt=line_split[6]
			ensg=line_split[50]
			TSS=int(line_split[44])-int(pos)
			TES=int(line_split[45])-int(pos)
			genetype=line_split[54]
			gnomad_split=(line_split[42]).split(';')
			vep_split=(line_split[42]).split('|')
			gene_start=int(line_split[44])
			gene_end=int(line_split[45])

			#this order matters -- will override if in both gene start and gene end window to record as gene start
			#again will override if actually within gene instead
			if(abs(int(pos)-gene_end) <= genewindow):
				var_location="gene_end"
			if(abs(int(pos)-gene_start) <= genewindow):
				var_location="gene_start"
			#if within gene
			if((int(pos)>gene_start) and (int(pos)<gene_end)):
				#get if exon or intron
				#while loop for speed.assuming sorting of both files.
				while(gtf_start<gene_start):
					gtf_line=(gtf.readline()).split("\t")
					gtf_start=int(gtf_line[3])
				gtf_end=int(gtf_line[4])
				if((int(pos)>gtf_start) and (int(pos)<gtf_end)):
					var_location="exon"
				else:
					var_location="intron"
			if((abs(int(pos)-gene_start) > genewindow) and (abs(int(pos)-gene_end) > genewindow)):
				continue

			###get vep anno. it repeats in lengths of ten. loop through that.
			num_vep_annos=int((len(vep_split)-1)/10)
			vepvar_set=set()
			snpeff_set=set()
			for i in range(0,num_vep_annos):
				this_vepvar=vep_split[i*10+1]
				this_vepsnpeff=vep_split[i*10+2]
				this_vepvar_split=this_vepvar.split("&")
				for myvepvar in this_vepvar_split:
					vepvar_set.add(myvepvar)
				snpeff_set.add(this_vepsnpeff)
			vepvar=','.join(str(s) for s in vepvar_set)
			vepsnpeff=','.join(str(s) for s in snpeff_set)

			print("passed: " +var_location)
			cadd_raw=line_split[len(line_split)-2]
			cadd_phred=(line_split[len(line_split)-1]).strip()
			##print(','.join([chr,pos,gtex_maf,ref,alt,ensg,cadd_raw,cadd_phred,genetype,str(TSS),str(TES),line_split[44],line_split[45]]))
			if gnomad_split[0] == "NO_MATCH":
				# print("NO MATCH: ",ref,",",alt,",",gtex_maf)
				#this prepares for cases where the reference allelse is actually teh minor allele
				if ref == alt:
					gnomad_maf_both=float(gtex_maf)
				#if it's been seen in gtex more than once, it shouldnt have MAF of zero					elif float(gtex_maf) > 0.001:
					gnomad_maf_both=float(gtex_maf)
				else:
					gnomad_maf_both=float(0)		
			else:
				vartype=(([col for col in gnomad_split if variant_type.match(col)])[0].split("="))[1]
				if(vartype!="snv"):
					gnomad_maf_both=gtex_maf
				#if not empty, make sure AF_nfe is actually found
				found_col=[col for col in gnomad_split if af_both.match(col)]
				if found_col:
					gnomad_maf_both=float((found_col[0].split("="))[1])
				else:
					gnomad_maf_both=gtex_maf
			#if gtex maf is large than 1% difference, use gtex maf
			if(float(gtex_maf)-float(gnomad_maf_both)>0.01):
				use_maf=gtex_maf
			else:
				use_maf=str(gnomad_maf_both)
			#if gnomad maf is zero, use gtex maf
			#this is bc we will be filtering for being in 2 people, so def not a real maf of zero...
			if(float(gnomad_maf_both)==0):
				use_maf=gtex_maf
			#print("\t".join([chr,pos,ensg, genetype,ind,sex,gtex_maf,str(gnomad_maf_both)]))
			myline=[chr,pos,pos,gtex_maf,str(gnomad_maf_both),use_maf,ind,"SNPs",ensg, genetype,sex,var_location,str(cadd_raw),str(cadd_phred),str(geno),str(TSS),str(TES),vepvar,vepsnpeff]
			# print(myline)
			#must be seen in both
			#length of dic will be two if has male and female!
			if (len(seen_dic[pos])<2):
				outf_sex_all.write(("\t".join(myline)+"\n").encode('utf-8'))
				continue
			else:
				outf_sex_inboth.write(("\t".join(myline)+"\n").encode('utf-8'))
				outf_sex_all.write(("\t".join(myline)+"\n").encode('utf-8'))
outf_sex_all.close()
outf_sex_inboth.close()
