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

#parser = argparse.ArgumentParser(description='argparser')
#parser.add_argument('dir_read', type=str, help='directory to read file from')
#parser.add_argument('cutoff_mafdiff', type=float, help='cutoff [0-1] for % difference')
#parser.add_argument('out_mafdiff', type=str, help='Output file for MAF Mvs.F diff')
#parser.add_argument('out_m', type=str, help='outfile collaped m')
#parser.add_argument('out_f', type=str, help='outfile collaped f')
#parser.add_argument('out_b', type=str, help='outfile collaped both')
#parser.add_argument('sex_file', type=str, help='GTEX sample file that can get sex')
#args = parser.parse_args()

dir_read="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/bySiteAnnoX/GenesAnno"
out_mafdiff="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/maf_diff_x.tsv.gz"
out_maf_both="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_both_x.tsv.gz"
out_maf_m="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_m_x.tsv.gz"
out_maf_f="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_f_x.tsv.gz"
sex_file="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"

#compile at the beginning for speed
af_both=re.compile("AF_nfe=")
af_m=re.compile("AF_nfe_male=")
af_f=re.compile("AF_nfe_female=")


outwrite_mafdiff=gzip.open(out_mafdiff,"wb")
outwrite_maf_both=gzip.open(out_maf_both,"wb")
outwrite_maf_m=gzip.open(out_maf_m,"wb")
outwrite_maf_f=gzip.open(out_maf_f,"wb")

#get sex per individual
sex_convert_key={}
sex_convert_key["1"]="male"
sex_convert_key["2"]="female"
print(sex_convert_key)
sex_key={}
with open(sex_file,"r") as sex_f_read:
	next(sex_f_read)
	for sex_line in sex_f_read.readlines():
		sex_line_split=sex_line.split("\t")
		sex_key[sex_line_split[0]]=sex_convert_key[sex_line_split[2]]
print("key for sex")
print(sex_key)


colname_mafdiff=["chr","pos","ensg"," vartype","ind","ref","alt","varswitch","gtex_maf",
			"gnomad_maf_both","gnomad_maf_m","gnomad_maf_f","genetype","gnomad_maf_diff"]
colname_collapsed=["chr","pos","ensg","vartype","ind","ref","alt","varswitch","gtex_maf","gnomad_maf_both",
			"gnomad_maf_m","gnomad_maf_f","genetype","gnomad_maf_diff","num_rvs"]
outwrite_mafdiff.write(('\t'.join(map(str,colname_mafdiff))+"\n").encode())
outwrite_maf_both.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())
outwrite_maf_m.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())
outwrite_maf_f.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())


### loop through SNPS/INDELS/SV
for vartype in ["SNP","indel","SV"]:
	## loop through all inds in this vartype
	for f in glob.glob(dir_read+"/*"+vartype+"*"):
		#open file and read line by line
		f_split=f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=f_split_again[1]
		sex=sex_key[ind]
		print(sex)

		#new dictionary for individuals
		#keys are genes
		#value is minim MAF
		dic_m=dict()
		dic_f=dict()
		dic_both=dict()


		with open(f,"r") as f_read:
			for line in f_read.readlines():
				line_split=line.split("\t")
				chr=line_split[0] #chr
				pos=line_split[1] #pos
				gtex_maf=line_split[3] #GTEx MAF
				ref=line_split[5] #reference allele
				alt=line_split[6] #alternate allele
				ensg=(line_split[42]).split('\"')[1] # ENSG ID
				genetype=((line_split[46]).split('\"'))[1] #lincRNA, proteincoding, etc
				varswitch=line_split[53].strip() ## A->C as AC
				gnomad_anno=line_split[34] # GNOMAD annotation to further split
				gnomad_split=gnomad_anno.split(';')
				gnomad_maf_both=float((([col for col in gnomad_split if af_both.match(col)])[0].split("="))[1])
				gnomad_maf_m=float((([col for col in gnomad_split if af_m.match(col)])[0].split("="))[1])
				gnomad_maf_f=float((([col for col in gnomad_split if af_f.match(col)])[0].split("="))[1])
	
				if (gnomad_maf_m+gnomad_maf_f)/2 == 0:
					gnomad_maf_diff=0
				else:
					gnomad_maf_diff=(gnomad_maf_m-gnomad_maf_f)/((gnomad_maf_m+gnomad_maf_f)/2)
				if gnomad_maf_diff > 0.1:
					outwrite_mafdiff.write(("\t".join([chr,pos,ensg, vartype,ind,ref,alt,varswitch,
						gtex_maf, str(gnomad_maf_both),str(gnomad_maf_m),str(gnomad_maf_f),
						genetype,str(gnomad_maf_diff)])+"\n").encode())
				#print("\t".join([chr,pos,ensg, vartype,ind,ref,alt,varswitch,
				#		gtex_maf, str(gnomad_maf_both),str(gnomad_maf_m),str(gnomad_maf_f),
				#		genetype,str(gnomad_maf_diff)]))

				if gnomad_maf_both < 0.01:
					store_line=[chr,pos,ensg, vartype,ind,ref,alt,varswitch,
						gtex_maf, gnomad_maf_both,gnomad_maf_m,gnomad_maf_f,
						genetype,gnomad_maf_diff]
					if ensg in dic_both:
						dic_current_min_maf=(dic_both[ensg])[9]
						if gnomad_maf_both < dic_current_min_maf:
							this_count=dic_both[ensg][-1]+1
							store_line.append(this_count)
							dic_both[ensg]=store_line
						else:
							dic_both[ensg][-1]+=1
					else:
						store_line.append(1)
						dic_both[ensg]=store_line
			##writes the dictionary to a file
			for key_both in dic_both:			
				outwrite_maf_both.write((key_both+"\t"+'\t'.join(map(str,dic_both[key_both]))+'\n').encode())
		break
	break				

outwrite_mafdiff.close()
outwrite_maf_both.close()
outwrite_maf_m.close()
outwrite_maf_f.close()

