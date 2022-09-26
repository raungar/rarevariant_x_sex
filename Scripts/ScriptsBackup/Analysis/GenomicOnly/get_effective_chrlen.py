import argparse
import glob
import re
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--dir_gtf', type=str, help='directory to read prot/lnc only genes')
parser.add_argument('--chr_lens', type=str, help='file with chr lens')
parser.add_argument('--outfile', type=str, help='outfile with effective chr lens')
args = parser.parse_args()

dir_gtf_f=args.dir_gtf
chr_lens_f=args.chr_lens
outfile=args.outfile

#dir_gtf_f="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/GTF_lnc_protcod_only"
#chr_lens_f="/oak/stanford/groups/smontgom/raungar/Sex/Files/GRCh38.p13_assembly_report.txt"
#outfile="tmp"


#make a dictionary with the number of bp in each chr
#the keys are the chr, the values are the chrlen
chr_len_dic={}
with open(chr_lens_f) as chr_len_read:
	next(chr_len_read) #skip header
	for chrlen_line in chr_len_read.readlines():
		chrlen_linesplit=chrlen_line.split("\t")
		this_chrlen=int(chrlen_linesplit[8])
		this_chr="chr"+str(chrlen_linesplit[0])
		chr_len_dic[this_chr]=this_chrlen
#effective chr len dictionary
#keys are chr, values are effective len
eff_len_dic={}

#loop through each gtf (chromosome)
for f in glob.glob(dir_gtf_f+"/*gtf"):
	this_chr=((f.split("/"))[-1]).split('.')[0]
	print(this_chr)
	#make a dictionary for this chromosome
	#key is all base pairs in this chr
	#val is initialized as zero, then 1 if any variant overlaps
	vars_dic={}
	this_chr_len=int(chr_len_dic[this_chr])
	for i in range(0,this_chr_len): 
		vars_dic[i]=0 
	with open(f) as this_f:
		for this_line in this_f.readlines():
			line_split=this_line.split("\t")
			start=int(line_split[1])-10000
			end=int(line_split[2])+10000
			#adjust if genes are toward end
			if(start<0):
				start=0
			if(end>this_chr_len):
				end=this_chr_len
			for j in range(start,end):
				vars_dic[j]=1
		this_effective_size=sum(vars_dic.values())
		print(str(this_effective_size))
		eff_len_dic[this_chr]=this_effective_size

with open(outfile,"w") as out:
	for k, v in eff_len_dic.items():
		out.write(str(k)+"\t"+str(v)+"\n")
