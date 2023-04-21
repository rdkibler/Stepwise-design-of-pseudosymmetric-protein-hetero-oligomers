#!/home/rdkibler/.conda/envs/pyro/bin/python3

import json
import glob

hbnet_paths = glob.glob("/net/scratch/rdkibler/201224_rotor_hbnet/*/*/*.pdb")
hbnet_mapping = {path.split("/")[-1].split(".pdb")[0] + "_designed":path for path in hbnet_paths}

with open("hbnet_mapping.json",'w') as f:
	json.dump(hbnet_mapping,f)

