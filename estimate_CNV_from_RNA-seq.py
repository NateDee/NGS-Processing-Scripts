# Create CNV calculator program from RNA-Seq data using methods from:
# https://www.ncbi.nlm.nih.gov/pubmed/24925914

# Start with sorted list of RPKM data, using genes with RPKM > 1
# Sort must be chromosome 1-22, X, Y, from low to high base pairs

# For each gene, CNV call = Summation of 50 prior, GOI, and 50 post genes in matrix; with each individual gene diveded by 101.
# Should be a simple script to write
# For chroms 1 and Y, only go to boundaries of chromosoms, don't pretend to circularize...
# Header of file should be, ['chr', 'start', 'end', 'transcriptID', 'GeneID', 'data columns.....']

import pandas as pd
import numpy as np
import sys

# Data should always start in column 6 (index 5)
rpkm_tmp = pd.read_table(sys.argv[1], sep="\t", header=0)

# Write checks to make sure formatting is correct:
rpkm = rpkm_tmp.sort_values(['chr', 'start']) # Make sure it's sorted
rpkm = rpkm.reset_index(drop=True)
# Write function for a gene call
def cnvGene(i, rpkm_np_col):  # i is the number of the gene in the matrix/list; rpkm_np_col is a numpy array of all sample gene rpkm values
	if i < 50:
		gene_range = range(0,i+51) # i+51 won't be included when iter through range
	elif i + 51 >= len(rpkm_np_col):
		gene_range = range(i-50, len(rpkm_np_col))
	else:
		gene_range = range(i-50,i+51)
	sum_exp = 0
	for x in gene_range:  # Update this for less than 101 observations
		sum_exp += (rpkm_np_col[x]/101)
	return(sum_exp)

# Function for processing each column of data, return pandas dataframe of final
def cnvRpkm(rpkmdf):
	data_range = range(5, len(rpkmdf.columns))
	res_list = [list() for i in range(0,len(data_range))]
	col_num = 0
	for col in data_range:
		rpkm_np_col = rpkmdf.iloc[:,col].values
		for gene in range(0, len(rpkm_np_col)):
			res_list[col_num].append(cnvGene(gene, rpkm_np_col))
		col_num += 1
	df_res_tmp = pd.DataFrame(res_list).transpose()
	df_res = pd.concat([rpkmdf.iloc[:, 0:5], df_res_tmp], axis = 1)
	df_res.columns = rpkmdf.columns
	return(df_res)
	
cnv_results = cnvRpkm(rpkm)
cnv_results.to_csv(sys.argv[2], sep="\t", header=True, index = False)

# After output manually process normalization (subtract average of your reference from the log2(RPKM) CNVcounts given here.  Zero center.
