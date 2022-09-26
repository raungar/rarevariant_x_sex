import argparse
import gzip

parser = argparse.ArgumentParser(description='argparser')
parser.add_argument('--outlierfile', type=str, help='directory to read file from')
parser.add_argument('--opentargetfile', type=str, help='directory to read file from')
parser.add_argument('--outfile', type=str, help='directory to read file from')
args = parser.parse_args()

gene2z={} #key is gene name (which is an outlier), value is zscore
with gzip.open(args.outlierfile,"rb") as outlierfile:
    for line in outlierfile:
        linespl=line.decode('utf-8').strip().split("\t")
        #if header
        if(linespl[0].startswith("Ind")):
            continue
        #only adding outliers to the dictionary
        if(linespl[5] != "outlier"):
            continue
        gene2z[linespl[1].split(".")[0]]=linespl[4]

outfile=open(args.outfile,"w")

#loop through open target file
with open(args.opentargetfile) as targetfile:
    for l in targetfile:
        l_spl=l.strip().split('\t')
        #in gene is one of the outliers
        this_gene=l_spl[1]
        if this_gene in gene2z:
            full_line=l_spl+[gene2z[this_gene]]+['\n']
            outfile.write('\t'.join(full_line))

outfile.close()
