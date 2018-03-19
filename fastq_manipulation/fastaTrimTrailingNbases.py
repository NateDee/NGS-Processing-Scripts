#This script will take a fasta file as input and output a fasta file that keeps only a specified number of bases, in this example, the last 1500 bases of each record.  If less than 1500 bases, will keep everything.

from Bio import SeqIO
import sys

for record in SeqIO.parse(sys.argv[1], "fasta"):
    try:    
        print(record[-1500:].format("fasta"))
    except ValueError:
        pass 
