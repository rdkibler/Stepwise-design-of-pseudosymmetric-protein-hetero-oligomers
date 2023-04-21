#!/bin/python3


import numpy as np



#num_points = 20
#pert_high = 5.0
#pert_low = 0.1

#for i in np.linspace(pert_low, pert_high, num_points):
#	print(f"./run/run.py inputs/RC4_20_mini.pdb 4 {i} 17 true")

for i  in range(250):
	print(f"./run/run.py inputs/RC4_20_mini.pdb {i} 1.6 31 true") #this creates 40 outputs, so if I want 40k, I just have to run 1000 of these! #let's start with 10k




