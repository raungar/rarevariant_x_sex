#!/bin/python

import argparse, re, gzip
#from collections import defaultdict


parser=argparse.ArgumentParser()
parser.add_argument("--bed_file", type=str, help="the bed file", required=True)
parser.add_argument("--zscore_file_both", type=str, help="the zscore file for males and females, final file from preprocessing snakefile", required=True)
parser.add_argument("--zscore_file_sex", type=str, help="the zscore file for this sex, final file from preprocessing snakefile", required=True)
parser.add_argument("--outfile_both", type=str, help="outfile for both group containing bed info and zscores", required=True)
parser.add_argument("--outfile_sex", type=str, help="outfile for this sex containing bed info and zscores", required=True)
parser.add_argument("--ind", type=str, help="the name of an individual (gtex-id)", required=True)
parser.add_argument("--sex", type=str, help="the sex of the individual", required=True)
args = parser.parse_args()


#make dictionary where key is gene, value is another dictionary
#value nested dictionary key is tissue type, value is z-score for given ind
dic_zscore_sex={} #defaultdict(dict)
with gzip.open(args.zscore_file_sex, 'rb') as zscore_f:
	loop_n=0
	for zline in zscore_f:
		zline_split=zline.decode('utf-8').split("\t")
		#for first line, which col is the gtex sample?
		if loop_n == 0 :
			ind_col=zline_split.index(str(args.ind))
			loop_n=1+loop_n
			continue
		#zline_split=zline.split("\t")
		#key does not exist, then make this as a list rather then appending
		if not zline_split[1]  in dic_zscore_sex:
			#make nested dictionary, value is  tissue key is ind_col
			this_line_dic={zline_split[0]:zline_split[ind_col]}
			#set the value of this key to a list
			dic_zscore_sex[zline_split[1]]=[]
			#append this line dictionary to the list
			dic_zscore_sex[zline_split[1]].append(this_line_dic)
		#if key already exists
		else:
			#set value dictionary
			this_line_dic={zline_split[0]:zline_split[ind_col]}
			#append value dictionary to gene key
			dic_zscore_sex[zline_split[1]].append(this_line_dic)
			#dic_zscore_sex[zline_split[2]]=dic_zscore_sex[zline_split[2]].append(list({dic_zscore_sex[zline_split[1]]:dic_zscore_sex[zline_split[ind_col]]}))
		loop_n=1+loop_n
		#dic_zscore_sex[zline_split[1]][zline_split[2]]=zline_split[ind_col]

###repeat for both
dic_zscore_both={} #defaultdict(dict)
in_both_group=True
with gzip.open(args.zscore_file_sex, 'rb') as zscore_f:
	loop_n=0
	for zline in zscore_f:
		zline_split=zline.decode('utf-8').split("\t")
		#for first line, which col is the gtex sample?
		if loop_n == 0 :
			ind_col=zline_split.index(str(args.ind))
			if not ind_col:
				print("This individual is not in the combined sex group")
				in_both_group=False
				break
			loop_n=1+loop_n
			continue
		#zline_split=zline.split("\t")
		#key does not exist, then make this as a list rather then appending
		if not zline_split[1]  in dic_zscore_both:
			#make nested dictionary, value is  tissue key is ind_col
			this_line_dic={zline_split[0]:zline_split[ind_col]}
			#set the value of this key to a list
			dic_zscore_both[zline_split[1]]=[]
			#append this line dictionary to the list
			dic_zscore_both[zline_split[1]].append(this_line_dic)
		#if key already exists
		else:
			#set value dictionary
			this_line_dic={zline_split[0]:zline_split[ind_col]}
			#append value dictionary to gene key
			dic_zscore_both[zline_split[1]].append(this_line_dic)
			#dic_zscore_both[zline_split[2]]=dic_zscore_both[zline_split[2]].append(list({dic_zscore_both[zline_split[1]]:dic_zscore_both[zline_split[ind_col]]}))
		loop_n=1+loop_n
		#dic_zscore_both[zline_split[1]][zline_split[2]]=zline_split[ind_col]
#print(dic_zscore)

len_tissues=len(dic_zscore_both["ENSG00000182378.13"])  #.values()) #[0]


outf_sex=gzip.open(args.outfile_sex,"wb")
if(in_both_group):
	outf_both=gzip.open(args.outfile_both,"wb")

with open(args.bed_file) as bed_f:
	loop_n=0;
	for bline in bed_f:
		bline_split=re.split("\t|;", bline)
		#print(bline_split)
		#chr/pos1/pos2/var/gtex_maf/gnomad_maf_all/gnomad_maf_f
		#print(len(bline_split))
		to_print_list=[bline_split[0],bline_split[1],bline_split[2],
				bline_split[len(bline_split)-1].strip(), bline_split[3],
				bline_split[36].split("=")[1], bline_split[108].split("=")[1], bline_split[124].split("=")[1]]
		#print(to_print_list)
		gene_index=[i for i,bline_split in enumerate(bline_split) if "gene_id" in bline_split ]#bline_split.index("gene_id") #196
		this_gene=bline_split[gene_index[0]].split("\"")[1]
		if this_gene in dic_zscore_sex.keys():
			tissue_dic_list=dic_zscore_sex[this_gene]
			#print(tissue_dic_list)
			tissue_list=[list(this_dic.values())[0] for this_dic in tissue_dic_list]
			#print(tissue_list)
		else:
			tissue_list=["NA"]*len_tissues
		
		descrip_list=[args.ind,args.sex,this_gene] #+tissue_list
		final_print_list=to_print_list+descrip_list+tissue_list
		outf_sex.write(('\t'.join(final_print_list)+"\n").encode('utf-8'))

		if in_both_group:
			if this_gene in dic_zscore_both.keys():
				tissue_dic_list=dic_zscore_both[this_gene]
				tissue_list_both=[list(this_dic.values())[0] for this_dic in tissue_dic_list]
			else:
	                     	tissue_list_both=["NA"]*len_tissues			

			final_print_list_both=to_print_list+descrip_list+tissue_list_both
			outf_both.write(('\t'.join(final_print_list_both)+"\n").encode('utf-8'))
		#print(bline_split[197]) #.split("=")[1])
#		#print ind1/chr1/pos/var1/var2/AC/sex/maf_all/maf_m/maf_f/gene/dic_zscore[gene].values()
#		#repeat for dic_zscore_difsex
#		
		
outf_sex.close()
outf_both.close()

