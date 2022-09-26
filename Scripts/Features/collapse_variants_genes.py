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
parser.add_argument('--dir_read', type=str, help='directory to read file from')
parser.add_argument('--cutoff_mafdiff', type=str, help='cutoff [0-1] for % difference')
parser.add_argument('--maf_filter',default=0.01, type=str, help='cutoff [0-1] for MAF of what is rare ')
parser.add_argument('--out_mafdiff', type=str, help='Output file for MAF Mvs.F diff')
parser.add_argument('--out_m', type=str, help='outfile collaped m')
parser.add_argument('--out_f', type=str, help='outfile collaped f')
parser.add_argument('--out_b', type=str, help='outfile collaped both')
parser.add_argument('--sex_file', type=str, help='GTEX sample file that can get sex')
args = parser.parse_args()

dir_read=args.dir_read
cutoff_mafdiff=float(args.cutoff_mafdiff)
min_maf=float(args.maf_filter)
out_mafdiff=args.out_mafdiff
out_maf_both=args.out_b
out_maf_m=args.out_m
out_maf_f=args.out_f
sex_file=args.sex_file

print("arguments parsed.")
#dir_read="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/bySiteAnnoX/GenesAnno"
#out_mafdiff="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/maf_diff_x.tsv.gz"
#out_maf_both="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_both_x.tsv.gz"
#out_maf_m="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_m_x.tsv.gz"
#out_maf_f="/oak/stanford/groups/smontgom/raungar/Sex/Output/features_v8/collapsed_maf_f_x.tsv.gz"
#sex_file="/oak/stanford/groups/smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS_v2_downloaded_april2020.txt"
#min_maf=0.01 #### maf to filter on 



#compile at the beginning for speed
af_both=re.compile("AF_nfe=")
af_m=re.compile("AF_nfe_male=")
af_f=re.compile("AF_nfe_female=")

##open files to write to
outwrite_mafdiff=gzip.open(out_mafdiff,"wb")
outwrite_maf_both=gzip.open(out_maf_both,"wb")
outwrite_maf_m=gzip.open(out_maf_m,"wb")
outwrite_maf_f=gzip.open(out_maf_f,"wb")

#get sex per individual
sex_convert_key={}
sex_convert_key["1"]="male"
sex_convert_key["2"]="female"
sex_key={}

#get sex of individual
with open(sex_file,"r") as sex_f_read:
	next(sex_f_read)
	for sex_line in sex_f_read.readlines():
		sex_line_split=sex_line.split("\t")
		sex_key[sex_line_split[0]]=sex_convert_key[sex_line_split[2]]

#write headers to files
colname_mafdiff=["chr","pos","ensg"," vartype","ind","sex","ref","alt","varswitch","gtex_maf",
			"gnomad_maf_both","gnomad_maf_m","gnomad_maf_f","genetype","gnomad_maf_diff"]
colname_collapsed=["chr","pos","ensg","vartype","ind","sex","ref","alt","varswitch","gtex_maf","gnomad_maf_both",
			"gnomad_maf_m","gnomad_maf_f","genetype","gnomad_maf_diff","num_rvs"]
outwrite_mafdiff.write(('\t'.join(map(str,colname_mafdiff))+"\n").encode())
outwrite_maf_both.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())
outwrite_maf_m.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())
outwrite_maf_f.write(('\t'.join(map(str,colname_collapsed))+"\n").encode())


### loop through SNPS/INDELS/SV
for vartype in ["SNP","indel","SV"]:
	## loop through all inds in this vartype
	for f in glob.glob(dir_read+"/*"+vartype+"*reduced.bed.gz"):
		print(f)
		#open file and read line by line
		#get individual ID and sex of this individual
		f_split=f.split("/")[-1]
		f_split_again=f_split.split("_")
		ind=next(filter(lambda x: re.search('GTEX',x),f_split_again))
		sex=sex_key[ind]
		## choose proper file to write to for sex specific
		if sex == "male":
			outwrite_sex=outwrite_maf_m			
		elif sex == "female":
			outwrite_sex=outwrite_maf_f
		else:
			print("ERROR: INVALID SEX")

		#new dictionary for individuals
		#keys are genes
		#value is minim MAF
		dic_sex=dict() # sex specific dictionary
		dic_both=dict()

		#actually open file, read line by line
		print("continue...")
		with gzip.open(f,"r") as f_read:
			for line in f_read.readlines():
				line_split=(line.decode('utf-8')).split("\t")
				#get interesting columns
				chr=line_split[0] #chr
				pos=line_split[1] #pos
				gtex_maf=line_split[3] #GTEx MAF
				geno=line_split[4]
				if(geno==0):
					continue
				ref=line_split[5] #reference allele
				alt=line_split[6] #alternate allele
				if(line_split[42] == "NA"):
					next
				ensg=(line_split[42]) #.split('\"')[1] # ENSG ID
				#genetype=((line_split[54]).split('\"'))[1] #lincRNA, proteincoding, etc
				genetype=((line_split[46])) #.split('\"'))[1] #lincRNA, proteincoding, etc
				#varswitch=line_split[61].strip() ## A->C as AC
				varswitch=line_split[53].strip() ## A->C as AC
				#gnomad_anno=line_split[42] # GNOMAD annotation to further split
				gnomad_anno=line_split[34] # GNOMAD annotation to further split
				gnomad_split=gnomad_anno.split(';')
				if gnomad_split[0] == "NO_MATCH":
					print("NO MATCH: ",ref,",",alt,",",gtex_maf)
					#this prepares for cases where the reference allelse is actually teh minor allele
					if ref == alt:
						gnomad_maf_both=float(gtex_maf)
						gnomad_maf_m=float(gtex_maf)
						gnomad_maf_f=float(gtex_maf)
					#if it's been seen in gtex more than once, it shouldnt have MAF of zero
					elif float(gtex_maf) > 0.001:
						gnomad_maf_both=float(gtex_maf)
						gnomad_maf_m=float(gtex_maf)
						gnomad_maf_f=float(gtex_maf)
					else:
						gnomad_maf_both=float(0)		
						gnomad_maf_m=float(0)		
						gnomad_maf_f=float(0)		
					print("both: ",gnomad_maf_both," , m: ", gnomad_maf_m, " , f: ", gnomad_maf_f)
				else:
					gnomad_maf_both=float((([col for col in gnomad_split if af_both.match(col)])[0].split("="))[1])
					gnomad_maf_m=float((([col for col in gnomad_split if af_m.match(col)])[0].split("="))[1])
					gnomad_maf_f=float((([col for col in gnomad_split if af_f.match(col)])[0].split("="))[1])
					if (float(gtex_maf)>0.01) and (gnomad_maf_both < 0.01): 
						gnomad_maf_both=float(gtex_maf)
					if (float(gtex_maf)>0.01) and (gnomad_maf_m < 0.01): 
						gnomad_maf_m=float(gtex_maf)
					if (float(gtex_maf)>0.01) and (gnomad_maf_f < 0.01): 
						gnomad_maf_f=float(gtex_maf)
				#get difference between MAF m vs F (controlling for divide by zero situation)
				if (gnomad_maf_m+gnomad_maf_f)/2 == 0:
					gnomad_maf_diff=0
				else:
					gnomad_maf_diff=(gnomad_maf_m-gnomad_maf_f)/((gnomad_maf_m+gnomad_maf_f)/2)
				#only write to file those with a sex MAF differences of 10%
				if gnomad_maf_diff > cutoff_mafdiff:
					outwrite_mafdiff.write(("\t".join([chr,pos,ensg, vartype,ind,sex,ref,alt,varswitch,
						gtex_maf, str(gnomad_maf_both),str(gnomad_maf_m),str(gnomad_maf_f),
						genetype,str(gnomad_maf_diff)])+"\n").encode())
				#print("\t".join([chr,pos,ensg, vartype,ind,sex,ref,alt,varswitch,
				#		gtex_maf, str(gnomad_maf_both),str(gnomad_maf_m),str(gnomad_maf_f),
				#		genetype,str(gnomad_maf_diff)]))

				#only store rare variants at a level of less than 0.01
				#this is the both case
				#if gnomad_maf_both < min_maf:
				if 1==1:
					store_line=[chr,pos,ensg, vartype,ind,sex,ref,alt,varswitch,
						gtex_maf, gnomad_maf_both,gnomad_maf_m,gnomad_maf_f,
						genetype,gnomad_maf_diff]
					#if there is already a RV recorded for this gene
					if ensg in dic_both:
						#get current min MAF
						dic_current_min_maf=float((dic_both[ensg])[9])
						#if this is more rare, store this instead
						if gnomad_maf_both < dic_current_min_maf:
							this_count=dic_both[ensg][-1]+1
							store_line.append(this_count)
							dic_both[ensg]=store_line
						else:
							# if this is not more rare, just inc the count of numRVs for this gene
							dic_both[ensg][-1]+=1
					#if there is not a RV recorded for this gene, store this one
					else:
						store_line.append(1)
						dic_both[ensg]=store_line
				#now do the same as above, but for the sex of this ind
				if sex == "male":
					gnomad_sex=gnomad_maf_m
				elif sex == "female":
					gnomad_sex=gnomad_maf_f
				else:
					print("ERROR: INVALID SEX")
				if gnomad_sex < min_maf:
					store_line=[chr,pos,ensg, vartype,ind,sex,ref,alt,varswitch,
						gtex_maf, gnomad_maf_both,gnomad_maf_m,gnomad_maf_f,
						genetype,gnomad_maf_diff]
					if ensg in dic_sex:
						dic_current_min_maf=float((dic_sex[ensg])[9])
						if gnomad_maf_both < dic_current_min_maf:
							this_count=dic_sex[ensg][-1]+1
							store_line.append(this_count)
							dic_sex[ensg]=store_line
						else:
							dic_sex[ensg][-1]+=1
					else:
						store_line.append(1)
						dic_sex[ensg]=store_line
			##writes the dictionary to a file
			for key_both in dic_both:			
				#outwrite_maf_both.write((key_both+"\t"+'\t'.join(map(str,dic_both[key_both]))+'\n').encode())
				outwrite_maf_both.write(('\t'.join(map(str,dic_both[key_both]))+'\n').encode())
			for key_both in dic_sex:			
				#outwrite_sex.write((key_both+"\t"+'\t'.join(map(str,dic_sex[key_both]))+'\n').encode())
				outwrite_sex.write(('\t'.join(map(str,dic_sex[key_both]))+'\n').encode())

outwrite_mafdiff.close()
outwrite_maf_both.close()
outwrite_maf_m.close()
outwrite_maf_f.close()
