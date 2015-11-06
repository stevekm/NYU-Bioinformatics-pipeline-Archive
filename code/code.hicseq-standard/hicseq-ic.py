#!/usr/bin/python

import sys
import numpy as np 
from mirnylib.numutils import completeIC

inp_file = sys.argv[1]
out_file = sys.argv[2]

inp_matrix = np.loadtxt(inp_file) 
out_matrix = completeIC(inp_matrix) 
np.savetxt(out_file,out_matrix,fmt='%.3f')

