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
		var_location=line_split[11]
		cadd_raw=line_split[12]
		cadd_phred=line_split[13].strip()
		geno=line_split[14].strip()
		if(cadd_phred=="NA"):
			cadd_phred=0
		#if homo dom, dont consider
		if(int(geno) == 0):
			continue
		store_line=[chrom,start,end,ensg, vartype,ind,sex,var_location, maf_gtex, maf_gnomad,maf_use,genetype,cadd_raw,cadd_phred,geno]
		#if there is already a RV recorded for this gene
		if ind not in this_inds_dic:
			this_inds_dic[ind]=dict()
		if ensg in this_inds_dic[ind]:
			#get current min MAF
			dic_current_min_maf=float((this_inds_dic[ind][ensg])[10])
			dic_current_min_cadd_phred=(this_inds_dic[ind][ensg])[13]
			dic_current_varloc=(this_inds_dic[ind][ensg])[7]
			#store if not intron!
			print(ind+":"+ensg+"---this/min maf: ",maf_use,"/",dic_current_min_maf, "   and this/mincadd ",
				cadd_phred,"/",dic_current_min_cadd_phred+". this/minvarloc=" +var_location+"/"+dic_current_varloc)
			#	this_count=this_inds_dic[ind][ensg][-1]+1
			#	store_line.append(this_count)
			#	this_inds_dic[ind][ensg]=store_line
			#	continue
			if((dic_current_varloc!="exon")  & (var_location =="exon")):
				print("storing bc not dic_current_varloc was " +dic_current_varloc+ " but this location is "+var_location)
				print(str(maf_use) + "=new vs. " +str(1.5*float(maf_use)) + "=1.5new vs. " +str(dic_current_min_maf)+"=old")
				#this is funky, but if it's an exon, only store if lower minor allele freq wihtin 1.5x
				if((float(maf_use)<float(dic_current_min_maf)*1.5) & (float(cadd_phred)>=cadd_min)):
					this_count=this_inds_dic[ind][ensg][-1]+1
					store_line.append(this_count)
					this_inds_dic[ind][ensg]=store_line
					print("storing bc of maf and exon")
					continue
			#if(var_location !="exon"):
			#		this_count=this_inds_dic[ind][ensg][-1]+1
			#		store_line.append(this_count)
			#		this_inds_dic[ind][ensg]=store_line
			#		continue
					
					
			#check first that this line passes the min cadd threshold (since will be sent to common eventually anyway)
			#and make sure th
			if(float(cadd_phred)>=cadd_min):
				#if this is more rare, store this instead
				#if float(maf_use) < dic_current_min_maf:
				#if  more rare...
				#please if there is an option that is reasonable that is not an intron grab it
				if(var_location == "intron"):
					if(dic_current_varloc != "intron"):
						print("plz no introns")
						#not an intron, so let's see if it's within 1.5 of more rare
						if(float(maf_use)<float(dic_current_min_maf)*1.5):					
							print("ok good no introns")
							this_count=this_inds_dic[ind][ensg][-1]+1
							store_line.append(this_count)
							this_inds_dic[ind][ensg]=store_line
							continue
				if((float(maf_use))<(float(dic_current_min_maf))):
					print(" more rare")
					this_count=this_inds_dic[ind][ensg][-1]+1
					store_line.append(this_count)
					this_inds_dic[ind][ensg]=store_line
					continue
				#store highest cadd
				if (float(dic_current_min_cadd_phred)<cadd_min):
					print("my cadd was rare")
					this_count=this_inds_dic[ind][ensg][-1]+1
					store_line.append(this_count)
					this_inds_dic[ind][ensg]=store_line
					continue # next if statement doesn't even matter, since needs to be stored regardless
			else:
				# if this is not more rare, just inc the count of numRVs for this gene
				this_inds_dic[ind][ensg][-1]+=1
				print("just increase count")
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
