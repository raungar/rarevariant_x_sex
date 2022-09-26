#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import gzip
import numpy as np
import os
############## FUNCTIONS

def sample_dict_maker(sample_fn):
	'''
	Reads in sample annotation file and makes dict. Keys are sample IDs, and values are tissues.
	'''
	sample_fh = open(sample_fn, 'r')
	sample_dict = {}
	for line in sample_fh:
		sample, tissue = line.rstrip().split('\t')
		if tissue != '':
			sample_dict[sample] = tissue
	sample_fh.close()
	return sample_dict

def split_by_tissue(sample_dict, gtex_fn, outdir, end):
	'''
	Read in GTEx TPM or reads file and split by tissues.
	'''
	gtex = gzip.open(gtex_fn, 'r')
	for line in gtex:
		line = line.decode('utf8').rstrip().split('\t')
		if 'Name' in line:
			header=line
			break
		if 'gene_id' in line:
			header = line
			break

	tissue_dict = {}
	for sample in header[2:]:
		if sample in sample_dict:
			tissue = sample_dict[sample]
			if tissue not in tissue_dict:
				tissue_dict[tissue] = []
			tissue_dict[tissue].append(header.index(sample))
	else:
		print(sample, 'not in sample-to-tissue map. Skipping it.', file = sys.stderr)

	header = np.array(header)
	file_dict = {}
	for tissue in tissue_dict:
		file_dict[tissue] = open(outdir + '/' + tissue + end, 'w')
		out_names = header[tissue_dict[tissue]]
		for i in range(len(out_names)):
			out_names[i] = '-'.join(out_names[i].split('-')[0:2])
		print('\t'.join(['Gene'] + out_names.tolist()), file = file_dict[tissue])

	for line in gtex:
		line = np.array(line.decode('utf8').rstrip().split())
		for tissue in tissue_dict:
			print('\t'.join([line[0]] + line[tissue_dict[tissue]].tolist()), file = file_dict[tissue])

	for tissue in file_dict:
		file_dict[tissue].close()



############## MAIN
print("split_expr_by_tissues.py")

usage = 'Splits combined RPKM file from GTEx into expression matrices by tissue.'

parser = argparse.ArgumentParser(description = usage)
parser.add_argument('--gtex', required = True, help = 'GTEx combined TPM or read count file')
parser.add_argument('--out', required = True, help = 'directory to hold output matrices for each tissue')
parser.add_argument('--sample', required = True, help = 'sample annotation file')
parser.add_argument('--end', required = True, help = 'file ending')

args = parser.parse_args()

###RAU at out if does not exist
if not os.path.exists(args.out):
	os.makedirs(args.out)

## Make sample dict and get list of tissues
sample_dict = sample_dict_maker(args.sample)

print(sample_dict)

## Split GTEx combined file by tissue
split_by_tissue(sample_dict, args.gtex, args.out, args.end)
