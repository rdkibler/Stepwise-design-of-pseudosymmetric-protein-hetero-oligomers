#!/home/rdkibler/.conda/envs/pyro/bin/python3

#didn't use silent files this time

#import sys
#sys.path.append("/home/rdkibler/scripts/silent_tools")
#import silent_tools

#silent_file = "/net/scratch/rdkibler/200629_churro_24x_normal_modes_perturb/churro_24x_normal_modes_perturbed_ASUs.silent"
#silent_index = silent_tools.get_silent_index(silent_file)

#decoys = silent_index['tags']


#for i,decoy in enumerate(decoys):
#	exe = f"./run/run.py -silent {silent_file} {decoy} -metavar production_varset.json -output_subdir {int(i/4000)}"
#	print(exe)


import glob
for i,path in enumerate(glob.glob("/net/scratch/rdkibler/201223_rotor_nmp/*/*.pdb")):
	print(f"./run/run.py -pdb {path} -metavar production_varset.json -output_subdir {int(i/4000)} ")
