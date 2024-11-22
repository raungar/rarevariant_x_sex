
# Takes PEER-corrected expression data 
# and combines it into a flat file.
# The output file is formatted as follows:
# 
#	tissue gene ind1_expr ind2_expr ...
#
# Missing values are coded as NAs
# Gets tissues and individual IDs from file
print(0)
import os
import numpy as np
from operator import itemgetter
import argparse

print(1)


parser = argparse.ArgumentParser()
parser.add_argument("-d","--RAREDIR",help="RARE DIR PATH, no trailing /")
parser.add_argument("-t","--tissue_file",help="path for tissue file (gtex_2017-06-05_tissues_all_normalized...)")
parser.add_argument("-i","--ind_file",help="path for individual file (gtex_2017-06-05_individuals_all_normalized...)")
parser.add_argument("-o","--outpath",help="path for gathered expression, include gz ending")
#parser.add_argument("-x","--x_gtf_file",help="path for x gtf file")
#parser.add_argument("-x","--a_gtf_file",help="path for autosomal gtf file")
parser.add_argument("-g","--group",help="group: sex.both, m, f")
args = parser.parse_args()



dir = args.RAREDIR #os.environ["RAREDIR"]
#dir = dir + '/preprocessing_v8/'
#tissueNamesFile = dir + 'gtex_2017-06-05_tissues_all_normalized_samples.txt'
#individualsFile = dir + 'gtex_2017-06-05_individuals_all_normalized_samples.txt'
#outfile = dir + 'gtex_2017-06-05_normalized_expression.txt'
tissueNamesFile = args.tissue_file
individualsFile=args.ind_file
outfile=args.outpath
my_group=args.group
#a_gtf_file=args.a_gtf_file
#x_gtf_file=args.x_gtf_file

exprdir = dir + '/PEER_v8/'


print("GATHERING:"+my_group)


# read IDs into a list
# add 'GTEX-' as a prefix and sort by id
ind = open(individualsFile, 'r')
individuals = [i.strip() for i in ind.readlines()]
individuals.sort()
ind.close()
#print("INDIVIDUALS")
#print(individuals)

# read in tissues to process in a list
# sort it
tis = open(tissueNamesFile, 'r')
tissues = [t.strip() for t in tis.readlines()]
tissues= [t.split("\t")[0] for t in tissues]
#print(tissues)
##del tissues[0] #remvoed rau 09/08/2020
#print(tissues)
tissues.sort()
#tissues.remove("Tissue") 
tissues=np.unique(tissues)
tis.close()
print("TISSUES")
print(tissues)

print("OUTFILE")
print(outfile)
# prepare output file
out = open(outfile, 'w')
out.write('\t'.join(['Tissue','Gene'] + individuals) + '\n')
#.decode('utf-8'))
# process tissues one at a time 
print("BEGINNING TO PROCESS")

for tissue in tissues:

	filename = exprdir + tissue + '.' + my_group+ '.peer.v8ciseQTL.ztrans.txt'
	print(filename)

	# read in header and figure out which columns to keep
	#try:
	if(1==1):
		expr = open(filename, 'r')
		headerList = expr.readline().strip().split()
		subjectIDs = headerList[1:]

		if(len(subjectIDs)<50):
			print("not adding "+ tissue+" because only "+str(len(subjectIDs))+ " individuals")
			continue
		cols2keep = [(i+1,subject) for (i,subject) in enumerate(subjectIDs) if subject in individuals]
        	# sort by subject ID (because brain_cerebellum wasn't sorted for some reason)
		sorted2keep = sorted(cols2keep, key=itemgetter(1))
		cols2keep = [0] + [s2k[0] for s2k in sorted2keep]
		index2add = [1] + [i+2 for (i,subject) in enumerate(individuals) if subject in subjectIDs]
		ngenes = sum(1 for _ in expr)
		expr.close()

		#print("NGENES: " , ngenes)

		# make sure that individual IDs line up
		orig = [(['Id'] + individuals)[i-1] for i in index2add ]
		new = [headerList[i] for i in cols2keep]
		assert len(new) == sum([o==n for o,n in zip(orig,new)])

		# initialize numpy array to contain final data
		# fill with NAs
		#padded = np.full((ngenes,len(individuals)+2), 'NA', dtype='|S40')
		#padded = np.full((ngenes,len(individuals)+2), 'NA', dtype='<U13')
		##padded = np.full((ngenes,len(individuals)+2), 'NA', dtype='<U100')
		padded = np.full((ngenes,len(individuals)+2),'NA', dtype='<U100')
		# add column with tissue
		padded[:,0] = tissue

		#print(0)
		# read data into numpy array
		values = np.loadtxt(filename, dtype=str, skiprows=1)
		# remove unwanted columns
		values = values[:, cols2keep]
		# put the values into the right slots in the preallocated array
		padded[:,index2add] = values
		#print("padded shape",padded.shape)
		#print(1)
		#padded2=np.array([x.decode() for x in padded])
		#print(2)
		# print matrix to file
		np.savetxt(out, padded,encoding='utf-8', fmt='%s', delimiter='\t', newline='\n')
	#except:
	#	print(tissue)

out.close()
