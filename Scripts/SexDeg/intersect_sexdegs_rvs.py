import argparse
import gzip
import re

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--genes_sexdegs', type=str, help='genes + sex deg file')
parser.add_argument('--rvsites_file', type=str, help='collpsed rv file')
parser.add_argument('--minmaf', type=float, help='min minor allel frequency')
parser.add_argument('--maxmaf', type=float, help='min minor allel frequency')
parser.add_argument('--outfile', type=str, help='outfile')
parser.add_argument('--logfile', type=str, help='logfile for ununsed genes')
args = parser.parse_args()
sexdegs_file=args.genes_sexdegs
rv_file=args.rvsites_file
outfile=args.outfile
logfile=args.logfile
minmaf=args.minmaf
maxmaf=args.maxmaf
outfile_write=gzip.open(outfile,"wb")
logfile_write=open(logfile,"w")

print("will filter between " + str(minmaf) + " and " + str(maxmaf))


#rv_file="Output/features_v8/collapsed_maf_both_x.tsv.gz"
#sexdegs_file="Output/sexdeg_v8/all_genes_sexDEGs_BREAST.txt"

#read in all genes into a dictionary, default beta (sexdeg) is 0
genes_dic={}
chr_dic={}

j=1

with open (sexdegs_file, "r") as f_degs:
	for line in f_degs.readlines():
		line_spl=line.strip().split("\t")
		print(line_spl)
		chr=line_spl[0]
		ensg=line_spl[1]
		beta=line_spl[2]
		chr_dic[ensg]=chr
		genes_dic[ensg]=[beta,0,1,"NA"] #this will be beta, inds, min MAF, sex [na/m/f/b]

i=1
with gzip.open(rv_file,"rb") as f_read:
	f_read.readline()
	for l in f_read.readlines():
		l_spl=l.decode('utf-8').strip().split("\t")
		print(l_spl)
		ensg=l_spl[2]
		sex=l_spl[5]
		maf=l_spl[10]
		print("i have: " +ensg+","+sex+","+maf)
		try: 
			gene_match=genes_dic[ensg]
		except:
			print("Key error thrown: "+ensg+" -- RV recorded but gene not included")
			logfile_write.write(ensg+"\n")
			continue
		gene_match_beta=float(gene_match[0])
		gene_match_ind=int(gene_match[1])
		gene_match_mafmin=float(gene_match[2])
		gene_match_sex=gene_match[3]

		print("maf: "+str(maf)+" , minmaf: "+str(minmaf))
		#SET MAF TO 1 if not within "RARE" def
		if (str(maf) == "gnomad_maf_both"):
			#THIS IS A HEADER -- SKIP
			print("THIS IS GNOMAD_MAF_BOTH")
			continue
		if not (float(maf) >= float(minmaf) and float(maf)<float(maxmaf)):
			maf=1
		print("i have: " +ensg+","+sex+","+str(maf)+" --- then.... "+str(gene_match_beta)+","+str(gene_match_sex)+","+str(gene_match_mafmin))

		#if gene_match_sex is NA , then no rv has been found near it yet, so skip ahead	
		if (gene_match_sex != "NA"):
			if (gene_match_mafmin !="NA"):
				if(float(gene_match_mafmin)<float(maf)):
					maf=gene_match_mafmin
			if (gene_match_sex != sex):
				sex="both"
		ind=gene_match_ind+1
		genes_dic[ensg]=[gene_match_beta,ind,maf,sex]
		####print('\t'.join(map(str,genes_dic[ensg])))


#finally, print the dictionary such that
#key is ensg, val is beta value
#print so tab delimeted ensg/beta

#for (k,v in genes_dic.iteritems()):
#for k,v in genes_dic.items():
#	outfile_write('\t'.join([str(k),map(str,v)]))

####coutfile_write.write('\t'.join([str(k),'\t'.join(map(str,v)]))
outfile_write.write(('\t'.join(["chr","ensg","beta","num_inds","min_maf","sex"])+"\n").encode('utf-8'))
for k,v in genes_dic.items():
	try:
		my_str=chr_dic[k]+"\t"+k+"\t"+'\t'.join(map(str,v))+"\n"
		outfile_write.write(my_str.encode('utf-8'))
	except:
		logfile_write.write(k+"\n")
	#print('\t'.join([str(k),map(str,v)]))

outfile_write.close()
logfile_write.close()
