import argparse
import glob
import re
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--gtf', type=str, help='gtf for prot/lnc only genes')
parser.add_argument('--chr_lens', type=str, help='file with chr lens')
parser.add_argument('--outfile', type=str, help='outfile with effective chr lens')
args = parser.parse_args()

gtf_f=args.gtf
chr_lens_f=args.chr_lens
outfile=args.outfile

#dir_gtf_f="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/GTF_lnc_protcod_only"
#chr_lens_f="/oak/stanford/groups/smontgom/raungar/Sex/Files/GRCh38.p13_assembly_report.txt"
#outfile="tmp"

chr_len_dic={}
#make a dictionary with the number of bp in each chr
#the keys are the chr, the values are the chrlen
#effective chr len dictionary
#keys are chr, values are effective len
par1_s=10001; par1_e=2781479; par2_s=155701383; par2_e=156030895; nonpar_s=2781480; nonpar_e=155701382;
xar_s=2731479;xar_e=58555579
xcr1_s=62462543; xcr1_e=89140830;
xtr_s=89140830;xtr_e=93428068;
xcr2_s=93428068;xcr2_e=155701383
len_par1=par1_e-par1_s;len_nonpar=par2_s-par1_e; len_par2=par2_e-par2_s
len_xcr1=xcr1_e-xcr1_s; len_xcr2=xcr2_e-xcr2_s
len_xtr=xtr_e-xtr_s; len_xar=xar_e-xar_s; 
chr_len_dic["PAR1"]=len_par1
chr_len_dic["PAR2"]=len_par2
chr_len_dic["NONPAR"]=len_nonpar
chr_len_dic["XCR1"]=len_xcr1
chr_len_dic["XCR2"]=len_xcr2
chr_len_dic["XAR"]=len_xar
chr_len_dic["XTR"]=len_xtr


eff_len_dic={}

#loop through each gtf (chromosome)
for this_key in chr_len_dic.keys():
	print(this_key)
	#make a dictionary for this chromosome
	#key is all base pairs in this chr
	#val is initialized as zero, then 1 if any variant overlaps
	vars_dic={}
	this_chr_len=int(chr_len_dic[this_key])
	for i in range(0,this_chr_len): 
		vars_dic[i]=0 
	f=gtf_f #glob.glob(digtf_f+"/*X*gtf")
	#this_chr=((f.split("/"))[-1]).split('.')[0]
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
		eff_len_dic[this_key]=this_effective_size

with open(outfile,"w") as out:
	for k, v in eff_len_dic.items():
		out.write(str(k)+"\t"+str(v)+"\n")
