### Using a strand aware sorted (- strand has lowest bp exon first, + strand has highest bp exon first) GTF file:
### Output modified gtf file that only has exons that encompass a given length of the 3' end
### Creating example for 1500 bp transcripts wanted

#Use pandas, then loop through dataframe
import numpy as np
import pandas as pd
import sys

#Read in gtf file to pandas dataframe

tmpdf = pd.read_table(sys.argv[1], sep="\t", header=None, names = ["chr", "genome", "type", "starts", "ends", "unused", "strand", "unused2", "gene_id"])



#Loop through pandas data frame
def trimGtf(df):
	results_df = pd.DataFrame(columns = ["chr", "genome", "type", "starts", "ends", "unused", "strand", "unused2", "gene_id"])
	exon_total = 0
	gene_done = "no"
	gene_unique = df['gene_id'].nunique()
	gene_count = 0
	for (i, row) in df.iterrows():
		#Columns are chr = 0, genome = 1, type = 2, starts = 3, ends = 4, unused = 5, strand = 6, unused2 =7, gene_id = 8
		gene = df.iloc[i,8]
		chrom = df.iloc[i,0]
		#Check for gene name in results_df
		#Check for gene
		if results_df['gene_id'].str.contains(gene).any():
			x_df = results_df[results_df['gene_id'].str.match(gene)].copy() # Get subset dataframe of results_df that match gene name
			x_df['difference'] = x_df['ends'] - x_df['starts'] # Calculate difference for all rows
			if x_df.loc[x_df['chr'] == chrom, 'difference'].sum() >= 1500: # If difference that matches gene name (and chromosome) is 1500, gene is done
				gene_done = "yes"
				exon_total = 1500
			elif x_df.loc[x_df['chr'] == chrom, 'difference'].sum() < 1500: # If difference is less than 1500, update exon_total
				gene_done = "no"
				exon_total = x_df.loc[x_df['chr'] == chrom, 'difference'].sum() # If chrom doesn't match exon_total will be 0 not NA so all good for that case
		else:
			gene_done = "no" # If not found in results, exon_total should be 0
			exon_total = 0
			gene_count += 1
			print(gene_count, "/", gene_unique, " genes done.", sep="")
		
		# If gene is not present, run the trimming
		if gene_done == "no":
			#Different if's for + / -
			if df.iloc[i,6] == "+":
				if (df.iloc[i,4] - df.iloc[i,3]) == 1500: #If exon length is 1500 add to dataframe
					results_df = results_df.append(df.iloc[i,])
					exon_total = 1500
				elif (df.iloc[i,4] - df.iloc[i,3]) < 1500: # If exon length is less than 1500, find bp_needed left
					bp_needed = 1500 - exon_total
					if (df.iloc[i,4] - df.iloc[i,3]) >= bp_needed: # If exon length > bp needed, add the slice of exon to results
						rowtmp = df.iloc[i,].copy()
						rowtmp.loc['starts'] = rowtmp.loc['ends'] - bp_needed
						results_df = results_df.append(rowtmp)
						exon_total = 1500
					elif (df.iloc[i,4] - df.iloc[i,3]) < bp_needed: # If exon lengh < bp needed, add the slice, update exon_total
						results_df = results_df.append(df.iloc[i,])
						# exon_total += df.iloc[i,4] - df.iloc[i,3] # Don't need this as exon total is counted each time
					
				elif (df.iloc[i,4] - df.iloc[i,3]) > 1500: # If exon length is greater than 1500, add just 1500 of that exon to dataframe
					rowtmp = df.iloc[i,].copy()
					rowtmp.loc['starts'] = rowtmp.loc['ends'] - 1500
					results_df = results_df.append(rowtmp)
					exon_total = 1500
				
			# If - strand	
			elif df.iloc[i,6] == "-":
				if (df.iloc[i,4] - df.iloc[i,3]) == 1500: #If exon length is 1500 add to dataframe
					results_df = results_df.append(df.iloc[i,])
					exon_total = 1500
				elif (df.iloc[i,4] - df.iloc[i,3]) < 1500: # If exon length is less than 1500, find bp_needed left
					bp_needed = 1500 - exon_total
					if (df.iloc[i,4] - df.iloc[i,3]) >= bp_needed: # If exon length > bp needed, add the slice of exon to results
						rowtmp = df.iloc[i,].copy()
						rowtmp.loc['ends'] = rowtmp.loc['starts'] + bp_needed
						results_df = results_df.append(rowtmp)
						exon_total = 1500
					elif (df.iloc[i,4] - df.iloc[i,3]) < bp_needed: # If exon lengh < bp needed, add the slice, update exon_total
						results_df = results_df.append(df.iloc[i,])
						# exon_total += df.iloc[i,4] - df.iloc[i,3]
					
				elif (df.iloc[i,4] - df.iloc[i,3]) > 1500: # If exon length is greater than 1500, add just 1500 of that exon to dataframe
					rowtmp = df.iloc[i,].copy()
					rowtmp.loc['ends'] = rowtmp.loc['starts'] + 1500
					results_df = results_df.append(rowtmp)
					exon_total = 1500
	return(results_df)	

results = trimGtf(tmpdf)
file_name = sys.argv[1] + "_1500bp_3prime.gtf"	
results.to_csv(file_name, sep="\t", header=False, index=False)
