#Simple script to count reads of fastq file provided as input.  Doing this with grep and command line is faster, btw.

from Bio import SeqIO
import sys

count = 0
for record in SeqIO.parse(sys.argv[1], "fastq"):
	count = count + 1

print("There are " + str(count) + " reads in file " + sys.argv[1])
