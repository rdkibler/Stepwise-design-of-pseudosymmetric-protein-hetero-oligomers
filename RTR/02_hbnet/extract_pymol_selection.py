#!/usr/bin/python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("pdb")

args = parser.parse_args()

with open(args.pdb,'r') as f:
	for line in f:
		if "pymol_selection" in line:
			sl = line.strip().split()
			name = "_".join(sl[0].split("_")[:-2])
			#print(name)
			#print(line)
			sl[2] = name + ","
			print(" ".join(sl[1:]))
