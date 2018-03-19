#Function to add 1 to all cells of text file
#Use this to add a constant before log-transforming data so that data is real numbers (i.e. no zeros being log transformed)
#Use this on a histgram (ghist) heatmap from HOMER, the result can then be passed to "cluster" program for clustering

import numpy as np
import pandas as pd
import sys

#Takes 1 argument, the filename of the heatmap.txt file - use a shell loop if you have multiple files
#Will overwrite that file

tmpfile = pd.read_table(sys.argv[1], sep="\t", index_col=0, header=0)
df = tmpfile + 1

df.to_csv(sys.argv[1], sep="\t", header=True, index=True)
