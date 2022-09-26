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
parser.add_argument('--combined_file', type=str, help='directory to read file from')
parser.add_argument('--cadd_min', type=str, help='directory to read file from')
parser.add_argument('--out', type=str, help='outfile collaped ')
args = parser.parse_args()

combined_file=args.combined_file
outfile=args.out
cadd_min=float(args.cadd_min)
print("arguments parsed.")



colnames_combined=["chr","start","end","maf_gtex","maf_gnomad","maf_use","ind","vartype","ensg","genetype","sex","cadd_raw","cadd_phred","geno"]

out=gzip.open(outfile,"wb")

#colname_collapsed=["chr","pos","ensg","vartype","ind","sex","gtex_maf","gnomad_maf","use_maf","genetype","num_rv"]
#out.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())

#key is ensg, value is number of rare variants (less than 0.01 ONLY) seen
seen_genes_rare=dict()
#keys is this ind, value is the other dict
this_inds_dic=dict()
count=0
with gzip.open(combined_file,"r") as f_read:
	for line in f_read.readlines():
		line_split=(line.decode('utf-8')).split("\t")
		#get interesting columns
		count+=1
		chrom=line_split[0]
		start=line_split[1]
		end=line_split[2]
		maf_gtex=line_split[3]
		maf_gnomad=line_split[4]
		maf_use=line_split[5]
		ind=line_split[6]
		vartype=line_split[7]
		ensg=line_split[8]
		genetype=line_split[9]
		sex=line_split[10]
		cadd_raw=line_split[11]
		cadd_phred=line_split[12].strip()
		geno=line_split[13].strip()
		if(cadd_phred=="NA"):
			cadd_phred=0
		#if homo dom, dont consider
		print(geno)
		print(int(geno) ==0)
		if(int(geno) == 0):
			continue
		store_line=[chrom,start,end,ensg, vartype,ind,sex, maf_gtex, maf_gnomad,maf_use,genetype,cadd_raw,cadd_phred,geno]
                out.write(('\t'.join(map(str,store_line))+'\n').encode())

out.close()
