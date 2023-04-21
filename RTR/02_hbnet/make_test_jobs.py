import numpy as np
from pathlib import Path
import json
import os

metavar_paths = []

base_varset = {"hb_threshold": -0.4, 
"onebody_hb_threshold": -0.2, 
"min_helices_contacted_by_network": 3, 
"min_network_size": 3, 
"max_unsat_Hpol": 0, 
"min_percent_hbond_capacity": 0.62, 
"min_intermolecular_hbonds": 2, 
"min_core_res": 3, 
"min_unique_networks": 1, 
"max_networks_per_pose": 3, 
"use_aa_dependent_weights": "true", 
"max_replicates_before_branch": 2, 
"max_replicates_before_unsat_check": 2, 
"only_symm_interfaces": "true", 
"at_least_one_net_w_aromatic_sc": "false", 
"seed_hbond_threshold": -1.2, 
"total_num_mc_runs": 10000, 
"store_subnetworks": "false", 
"interface_distance": 10.0, 
"core_cutoff": 4.2, 
"deep_core_cutoff": 5.2, 
"surface_cutoff": 1.8, 
"deep_surface_cutoff": 2.0}

Path("test_vars").mkdir(parents=True, exist_ok=True)
for hb_threshold in [-0.8, -0.4]:
	for onebody_hb_threshold in [-0.6, -0.2]:
		for max_unsat_Hpol in [0,1]:
			for min_percent_hbond_capacity in [0.62, 0.8]:
				for max_replicates_before_branch in [0,2]:
					for max_replicates_before_unsat_check in [0,2]:
						for total_num_mc_runs in [1000]:
							for interface_distance in [7.0, 10.0]:
								varset_name = f"varset_{hb_threshold}_{onebody_hb_threshold}_{max_unsat_Hpol}_{min_percent_hbond_capacity}_{max_replicates_before_branch}_{max_replicates_before_unsat_check}_{total_num_mc_runs}_{interface_distance}.vars"
								varset = dict(base_varset)
								varset["hb_threshold"] = hb_threshold
								varset["onebody_hb_threshold"] = onebody_hb_threshold
								varset["max_unsat_Hpol"] = max_unsat_Hpol
								varset["min_percent_hbond_capacity"] = min_percent_hbond_capacity
								varset["max_replicates_before_branch"] = max_replicates_before_branch
								varset["max_replicates_before_unsat_check"] = max_replicates_before_unsat_check
								varset["total_num_mc_runs"] = total_num_mc_runs
								varset["interface_distance"] = interface_distance
								with open("test_vars/"+varset_name,'w') as f:
									json.dump(varset,f)
								metavar_paths.append("test_vars/"+varset_name)


#just randomly picked decoys
decoy_list = [
"churro-24x_nm_02939_2.0_60",
"churro-24x_nm_00374_2.0_60",
"churro-24x_nm_01475_2.0_60",
"churro-24x_nm_01885_2.0_60",
"churro-24x_nm_00539_2.0_60",
"churro-24x_nm_00372_2.0_60",
"churro-24x_nm_02082_2.0_60",
"churro-24x_nm_00209_2.0_60",
"churro-24x_nm_01138_2.0_60",
"churro-24x_nm_00469_2.0_60",
]


for metavar_path in metavar_paths:
	for decoy in decoy_list:

		cmd = f"./run/run.py " \
		      f"-silent /net/scratch/rdkibler/200629_churro_24x_normal_modes_perturb/churro_24x_normal_modes_perturbed_ASUs.silent {decoy} " \
		      f"-debug " \
		      f"-metavar {metavar_path} " \
		      f"-output_suffix _{os.path.splitext(os.path.basename(metavar_path))[0]} " \

		print(cmd)
