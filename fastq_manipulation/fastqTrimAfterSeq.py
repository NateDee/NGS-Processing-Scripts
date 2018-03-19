#This script will take a fastq file as input and output a fastq file that is trimmed before a specified sequence defined by positionSTART variable. Can be changed to make positionEND to trim anything after a specified sequence.

from Bio import SeqIO
import sys

for record in SeqIO.parse(sys.argv[1], "fastq"):
    try:    
        positionSTART = str(record.seq).index("GGGGCCATGG")
        print(record[positionSTART:].format("fastq"))
    except ValueError:
        pass 
