print("EXCUSE ME")

import argparse
import gzip
import os
import re
import datetime
import glob
print("UMMM")
parser = argparse.ArgumentParser(description='argparser')
parser.add_argument("--combined_file", type=str, help="the bed file", required=True)
parser.add_argument('--outlierfile', type=str, help='directory to read file from')
parser.add_argument('--outfile_rawnum', type=str, help='directory to read file from')
parser.add_argument('--maf_min', type=str, help='directory to read file from')
parser.add_argument('--cadd_min', type=str, help='directory to read file from')
parser.add_argument('--maf_max', type=str, help='directory to read file from')
parser.add_argument('--outfile', type=str, help='directory to read file from')
args = parser.parse_args()
print("HELLLLO")
outliers=[]
nonoutliers=[]
i=0
j=0
print("here")
print(args.outlierfile)
#store which genes are outliers and which are not
print(datetime.datetime.now())
with gzip.open(args.outlierfile,"rb") as outlier_f:
	print("opened")
	for line in outlier_f:
		i=i+1
		if(i==1):
			#skip header
			continue
		line_split=line.decode('utf-8').split("\t")
		ind=line_split[0]
		gene=line_split[1]
		is_outlier=line_split[5]
		if(is_outlier=="outlier"):
			outliers.append({ind:gene})
			j=j+1
			print("outlier")
		else:
			nonoutliers.append({ind:gene})

		# if(i==20):
		# 	break
# print(outliers)
# print(nonoutliers)
print(j)
print("OK")     
print(len(outliers))  
print("appended")
print(datetime.datetime.now())

outliers_set=set()
nonoutliers_set=set()
outlier_dic={}
nonoutlier_dic={}
count=0
with gzip.open(args.combined_file,"r") as f_read:
	for line in f_read.readlines():
		line_split=(line.decode('utf-8')).split("\t")
		#get interesting columns
		count+=1		
		geno=line_split[14].strip() #if homo dom, dont consider
		if(int(geno) == 0):
			continue
		maf_use=line_split[5]
		#don't consider anything not in maf range
		if(not(maf_use>=args.maf_min and maf_use<args.maf_max)):
			continue
		cadd_phred=line_split[13].strip()
		if(cadd_phred=="NA"):
			cadd_phred=0
		if(float(cadd_phred)<float(args.cadd_min)):
			continue
		this_ind=line_split[6]
		# vartype=line_split[7]
		ensg=line_split[8]
		veptype=line_split[17]
		veptype_spl=veptype.split(",")
		# vepconsq=line_split[18].strip()
		if({this_ind:ensg} in outliers):
			outliers_set.add(this_ind+','+ensg)
			for this_vep in veptype_spl:
				if(this_vep not in nonoutlier_dic):
					outlier_dic[(this_ind,ensg,this_vep)]=1
		# elif({this_ind:ensg in nonoutliers}):
		else:
			nonoutliers_set.add(this_ind+','+ensg)
			for this_vep in veptype_spl:
				if(this_vep not in nonoutlier_dic):
					nonoutlier_dic[(this_ind,ensg,this_vep)]=1
		# else:
		# 	print("ERROR: not in outlier or nonoutlier list!!!")
		# if(len(nonoutliers_set)==5):
		# 	break
print(datetime.datetime.now())

# print(outliers not in )
print(len(nonoutliers_set))
# outlier_dic_prop={(k2):v/len(outliers_set) for k1,k2,k3 in k for k,v in outlier_dic.items()} #divide by num of gene x ind outliers 
gene_ind_pairs_nonoutlier=set()
nonoutlier_dic_vartypes={}
for k,v in nonoutlier_dic.items():
	gene_ind_pairs_nonoutlier.add((k[0],k[1]))
	if k[2] not in nonoutlier_dic_vartypes:
		nonoutlier_dic_vartypes[k[2]]=1
	else:
		nonoutlier_dic_vartypes[k[2]]+=1
print(datetime.datetime.now())

gene_ind_pairs_outlier=set()
outlier_dic_vartypes={}
for k,v in outlier_dic.items():
	gene_ind_pairs_outlier.add((k[0],k[1]))
	if k[2] not in outlier_dic_vartypes:
		outlier_dic_vartypes[k[2]]=1
	else:
		outlier_dic_vartypes[k[2]]+=1

print(datetime.datetime.now())
print("mathing")
nonoutlier_dic_prop={k:v/len(gene_ind_pairs_nonoutlier) for k,v in nonoutlier_dic_vartypes.items()}
outlier_dic_prop={k:v/len(gene_ind_pairs_outlier) for k,v in outlier_dic_vartypes.items()}

gene_name_outliers_found=[{x.split(',')[0]:x.split(',')[1]} for x in outliers_set]
outliers_norvs=[allout not in gene_name_outliers_found for allout in outliers]
outlier_dic_vartypes["no_variant"]=len([x for x in outliers_norvs if x])
outlier_dic_prop["no_variant"]=outlier_dic_vartypes["no_variant"]/len(outliers)

gene_name_nonoutliers_found=[{x.split(',')[0]:x.split(',')[1]} for x in nonoutliers_set]
nonoutliers_norvs=[allout not in gene_name_nonoutliers_found for allout in nonoutliers]
nonoutlier_dic_vartypes["no_variant"]=len([x for x in nonoutliers_norvs if x])
nonoutlier_dic_prop["no_variant"]=nonoutlier_dic_vartypes["no_variant"]/len(nonoutliers)
print(datetime.datetime.now())

print("write0")
with open(args.outfile_rawnum,"w") as outfile_rawnum:
	for k,v in outlier_dic_vartypes.items():
		outfile_rawnum.write('\t'.join([k,str(v),"outlier",'\n']))
	for k,v in nonoutlier_dic_vartypes.items():
		outfile_rawnum.write('\t'.join([k,str(v),"nonoutlier",'\n']))
print(datetime.datetime.now())
print("write")
with open(args.outfile,"w") as outfile:
	for k,v in outlier_dic_prop.items():
		outfile.write('\t'.join([k,str(v),"outlier",'\n']))
	for k,v in nonoutlier_dic_prop.items():
		outfile.write('\t'.join([k,str(v),"nonoutlier",'\n']))

print(datetime.datetime.now())
