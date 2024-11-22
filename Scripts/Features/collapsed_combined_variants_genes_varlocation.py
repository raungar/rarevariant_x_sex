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
parser.add_argument('--cadd_min', type=str, help='min cadd score')
parser.add_argument('--vartype', type=str, help='do this by each vartype')
parser.add_argument('--out', type=str, help='outfile collaped ')
args = parser.parse_args()

combined_file=args.combined_file
outfile=args.out
filtervar=args.vartype
cadd_min=float(args.cadd_min)
print("arguments parsed.")



#colnames_combined=["chr","start","end","maf_gtex","maf_gnomad","maf_use","ind","vartype","ensg","genetype","sex","var_location","cadd_raw","cadd_phred","geno"]

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
		var_location=line_split[11]
		#check if this matches the variant of interest
		if(var_location != filtervar):
			continue
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
		cadd_raw=line_split[12]
		cadd_phred=line_split[13].strip()
		geno=line_split[14].strip()
		veptype=line_split[17]
		vepconsq=line_split[18].strip()
		if(cadd_phred=="NA"):
			cadd_phred=0
		#if homo dom, dont consider
		if(int(geno) == 0):
			continue
		store_line=[chrom,start,end,ensg, vartype,ind,sex,var_location, maf_gtex, maf_gnomad,maf_use,genetype,cadd_raw,cadd_phred,geno,veptype,vepconsq]
		print(ensg)
		#if there is already a RV recorded for this gene
		if ind not in this_inds_dic:
			this_inds_dic[ind]=dict()
		if ensg in this_inds_dic[ind]:
			#get current min MAF
			dic_current_min_maf=float((this_inds_dic[ind][ensg])[10])
			dic_current_min_cadd_phred=(this_inds_dic[ind][ensg])[13]

			print("this/min maf: ",maf_use,"/",dic_current_min_maf, "   and this/mincadd ",cadd_phred,"/",dic_current_min_cadd_phred)
			#check first that this line passes the min cadd threshold (since will be sent to common eventually anyway)
			#and make sure th
			if(float(cadd_phred)>=cadd_min):
				#if this is more rare, store this instead
				if float(maf_use) < dic_current_min_maf:
					this_count=this_inds_dic[ind][ensg][-1]+1
					store_line.append(this_count)
					this_inds_dic[ind][ensg]=store_line
					continue
				#store highest cadd
				if (float(dic_current_min_cadd_phred)<cadd_min):
					this_count=this_inds_dic[ind][ensg][-1]+1
					store_line.append(this_count)
					this_inds_dic[ind][ensg]=store_line
					continue # next if statement doesn't even matter, since needs to be stored regardless
			else:
				# if this is not more rare, just inc the count of numRVs for this gene
				this_inds_dic[ind][ensg][-1]+=1
		#if there is not a RV recorded for this gene, store this one
		else:
			store_line.append(1)
			this_inds_dic[ind][ensg]=store_line

##writes the dictionary to a file
for this_ind in this_inds_dic:
	for this_gene in this_inds_dic[this_ind]:
		this_line=this_inds_dic[this_ind][this_gene]
		# print('\t'.join(map(str,this_line))+'\n')
		out.write(('\t'.join(map(str,this_line))+'\n').encode())
out.close()
