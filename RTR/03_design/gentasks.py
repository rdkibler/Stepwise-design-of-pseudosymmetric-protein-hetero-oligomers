#!/usr/bin/python3
import glob


#pdbs = glob.glob("/net/scratch/rdkibler/210121_rotor_rebuild_hbnet/*/*/*.pdb")

with open("../03_hbnet/save_passing.txt",'r') as f:
	pdbs = [line.strip() for line in f]

max_elements_per_dir = 2000
elements_in_dir = 0
dir_counter = 0

for pdb in pdbs:
	if elements_in_dir > max_elements_per_dir:
		elements_in_dir = 0
		dir_counter += 1
	elements_in_dir += 1
	print(f"./run/run.py {pdb} {dir_counter}")
