import glob


pdbs = glob.glob("../02_close_loop/output/*/*.pdb")

max_elements_per_dir = 2000
elements_in_dir = 0
dir_counter = 0

for pdb in pdbs:
	if elements_in_dir > max_elements_per_dir:
		elements_in_dir = 0
		dir_counter += 1
	elements_in_dir += 1
	print(f"./run/run.py {pdb} {dir_counter}")
