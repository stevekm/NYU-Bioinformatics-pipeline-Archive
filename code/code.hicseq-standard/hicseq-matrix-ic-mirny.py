#!/usr/bin/python

import sys
import numpy as np 
from mirnylib.numutils import completeIC

inp_file = sys.argv[1]
percent_cutoff = float(sys.argv[2])
out_file = sys.argv[3]

# read matrix
inp_matrix = np.loadtxt(inp_file) 

# filter matrix by row mean
m = np.mean(inp_matrix,axis=1)
cutoff = np.sort(m)[len(m)*percent_cutoff]
i = np.where(m<cutoff)
inp_matrix[i,:] = 0
inp_matrix[:,i] = 0

# IC
out_matrix = completeIC(inp_matrix) 

np.savetxt(out_file,out_matrix,fmt='%.3e')          # TODO: save in scientific notation

