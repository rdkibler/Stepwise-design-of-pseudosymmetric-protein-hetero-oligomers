#!/usr/bin/python3

import uuid
import sys
import subprocess
import os
import re
import shlex
import json
import time
import argparse

parser = argparse.ArgumentParser()

struct_control = parser.add_mutually_exclusive_group(required=True)
struct_control.add_argument("-pdb",metavar='PATH',type=str,help="path to input pdb")
struct_control.add_argument("-silent",metavar=('PATH','TAG'),nargs=2,help="path to silent file and a tag that exists in that file")

parser.add_argument("-metavars",type=str,default="",help="path to json file containing script_vars")
parser.add_argument("-output_suffix",type=str,default="")
parser.add_argument("-output_subdir",type=str,default=None)
parser.add_argument("-debug",action="store_true",default=False)

hbnet_control = parser.add_argument_group("HBNet controls")
hbnet_control.add_argument("-hb_threshold",type=float,default=-0.1,help="default: %(default)s") # %%hb_threshold%%
hbnet_control.add_argument("-onebody_hb_threshold",type=float,default=-0.4,help="default: %(default)s") # %%onebody_hb_threshold%%
hbnet_control.add_argument("-min_helices_contacted_by_network",type=int,default=2,help="default: %(default)s") # %%min_helices_contacted_by_network%%
hbnet_control.add_argument("-min_network_size",type=int,default=3,help="default: %(default)s") # %%min_network_size%%
hbnet_control.add_argument("-max_unsat_Hpol",type=int,default=0,help="default: %(default)s") # %%max_unsat_Hpol%%
hbnet_control.add_argument("-no_heavy_unsats_allowed",choices=['true','false'],default='false',help="default: %(default)s")
hbnet_control.add_argument("-min_percent_hbond_capacity",type=float,default=0.5,help="default: %(default)s") # %%min_percent_hbond_capacity%%
hbnet_control.add_argument("-min_intermolecular_hbonds",type=int,default=1,help="default: %(default)s") # %%min_intermolecular_hbonds%%
hbnet_control.add_argument("-min_core_res",type=int,default=2,help="default: %(default)s") # %%min_core_res%%
hbnet_control.add_argument("-min_unique_networks",type=int,default=1,help="default: %(default)s") # %%min_unique_networks%%
hbnet_control.add_argument("-max_networks_per_pose",type=int,default=3,help="default: %(default)s") # %%max_networks_per_pose%%
hbnet_control.add_argument("-use_aa_dependent_weights",choices=['true','false'],default='true',help="default: %(default)s") # %%use_aa_dependent_weights%%
hbnet_control.add_argument("-verbose",choices=['true','false'],default='true',help="default: %(default)s") # %%verbose%%
hbnet_control.add_argument("-max_replicates_before_branch",type=int,default=0,help="default: %(default)s") # %%max_replicates_before_branch%% ##this was 1 before
hbnet_control.add_argument("-max_replicates_before_unsat_check",type=int,default=0,help="default: %(default)s") # %%max_replicates_before_branch%% ##this was 1 before
hbnet_control.add_argument("-only_symm_interfaces",choices=['true','false'],default='false',help="default: %(default)s")
hbnet_control.add_argument("-at_least_one_net_w_aromatic_sc",choices=['true','false'],default='false',help="default: %(default)s")
hbnet_control.add_argument("-monte_carlo_seed_must_be_fully_buried",choices=['true','false'],default='false',help="default: %(default)s")

hbnet_control.add_argument("-only_start_at_interface_pairs",choices=['true','false'],default='true',help="Only start IG traversal with h-bonds that span across interface (different chains). default: %(default)s")

hbnet_control.add_argument("-seed_hbond_threshold",type=float,default=-0.4,help="default: %(default)s") # %%seed_hbond_threshold%%
hbnet_control.add_argument("-total_num_mc_runs",type=int,default=100000,help="default: %(default)s") # %%total_num_mc_runs%%
hbnet_control.add_argument("-store_subnetworks",choices=['true','false'],default='false',help="default: %(default)s") # %%store_subnetworks%% ##this was true before and I think that was a mistake
hbnet_control.add_argument("-interface_distance",type=float,default=8.0,help="default: %(default)s")

layer_control = parser.add_argument_group("layer controls")
layer_control.add_argument("-core_cutoff",type=float,default=4.2,help="default: %(default)s") # %%core_co%%
layer_control.add_argument("-deep_core_cutoff",type=float,default=5.2,help="default: %(default)s") # %%deep_core_cutoff%%
layer_control.add_argument("-surface_cutoff",type=float,default=1.8,help="default: %(default)s") # %%surf_co%%
layer_control.add_argument("-deep_surface_cutoff",type=float,default=2.0,help="default: %(default)s") # %%deep_surf_co%%

args = parser.parse_args()


args_as_dict = vars(args)
script_vars_args = set(["hb_threshold","onebody_hb_threshold","min_helices_contacted_by_network","min_network_size","max_unsat_Hpol","min_percent_hbond_capacity",
					"min_intermolecular_hbonds","min_core_res","min_unique_networks","max_networks_per_pose","use_aa_dependent_weights",
					"max_replicates_before_branch","max_replicates_before_unsat_check","seed_hbond_threshold","total_num_mc_runs","store_subnetworks","core_cutoff",
					"deep_core_cutoff","surface_cutoff","deep_surface_cutoff","interface_distance","only_symm_interfaces","at_least_one_net_w_aromatic_sc",
					"monte_carlo_seed_must_be_fully_buried","only_start_at_interface_pairs","no_heavy_unsats_allowed"])

script_vars_dict = {k:v for k,v in args_as_dict.items() if k in script_vars_args}

def format_script_vars(svd):
	return " ".join([f"{k}={str(v)}" for k,v in svd.items()])

#production = not args.debug





def run_command(command, outfile,to_file):
	if to_file:
		with open(outfile, 'w+') as out:
			try:
				dig = os.environ["SLURMD_NODENAME"]
				try:
					if os.environ["SLURM_ARRAY_TASK_ID"] == "":
						job = os.environ["SLURM_JOB_ID"]
					else:
						job = os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]
				except KeyError: 
					job = os.environ["SLURM_JOB_ID"]
			except KeyError:
				dig = "tack"
				job = "local"

			out.write("running on %s with jobid %s\n" % (dig, job))
			out.flush()

			start_time = time.time()
			result = subprocess.run(shlex.split(command), stdout=out, stderr=out)
			end_time = time.time()
			out.write(f"\nJOB COMPLETED IN {end_time - start_time}s\n")
			out.flush()
			#result = subprocess.run(shlex.split(command), stdout=subprocess.PIPE)
			#out.write(result.stdout.decode('utf-8'))
	else:
		start_time = time.time()
		subprocess.run(shlex.split(command))
		end_time = time.time()
		print(f"JOB COMPLETED IN {end_time - start_time}s")
		#with open(outfile, 'w+') as out:
		#	out.write(f"JOB COMPLETED IN {end_time - start_time}s\n")


# def run_command(command, outfile, to_file=False, to_screen=False):

# 	with open(outfile, 'w+') as out:
# 		try:
# 			dig = os.environ["SLURMD_NODENAME"]
# 			try:
# 				if os.environ["SLURM_ARRAY_TASK_ID"] == "":
# 					job = os.environ["SLURM_JOB_ID"]
# 				else:
# 					job = os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]
# 			except KeyError: 
# 				job = os.environ["SLURM_JOB_ID"]
# 		except KeyError:
# 			dig = "tack"
# 			job = "local"

# 		out.write("running on %s with jobid %s\n" % (dig, job))
# 		out.flush()

# 		start_time = time.time()


# 		#result = subprocess.run(shlex.split(command), stdout=out, stderr=out)
# 		proc = subprocess.Popen(shlex.split(command), 
# 								stdout=subprocess.PIPE, 
# 								stderr=subprocess.PIPE)


# 		while proc.poll() is None:
# 			line = proc.stderr.readline()
# 			if line:
# 				line = line.decode("utf-8")
# 				print("err: " + line,end="")
# 				out.write(line)
# 			line = proc.stdout.readline()
# 			if line:
# 				line = line.decode("utf-8")
# 				print("out: " + line,end="")
# 				out.write(line)

# 		end_time = time.time()
# 		out.write(f"\nJOB COMPLETED IN {end_time - start_time}s\n")
# 		out.flush()
# 		#result = subprocess.run(shlex.split(command), stdout=subprocess.PIPE)
# 		#out.write(result.stdout.decode('utf-8'))





if args.pdb:
	name = os.path.splitext(os.path.basename(args.pdb))[0]
else:
	name = args.silent[1]

name += args.output_suffix

outdir = "/net/scratch/rdkibler/201224_rotor_hbnet"

if args.output_subdir:
	outdir += "/" + args.output_subdir

outpath = os.path.abspath(outdir + "/%s" % name)

from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)


if args.metavars != "":
	#load these instead of the commandline ones
	with open(args.metavars,'r') as f:
		script_vars_dict = json.load(f)


with open(f"{outpath}/{name}_hb_metavars.txt",'w') as f:
	json.dump(script_vars_dict,f)


protocol="xml/hbnet.xml"

rosetta="/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccrelease"

if args.pdb:
	input_lines = f"-s {args.pdb} "
else:
	input_lines = f"-in:file:silent_struct_type binary -in:file:silent {args.silent[0]} -in:file:tags {args.silent[1]} "

extra_lines = ""
if args.debug:
	extra_lines = "-constant_seed " \
				  "-jran 987123 " \

	script_vars_dict['verbose'] = "true"
else:
	extra_lines = "-mute all "
	script_vars_dict['verbose'] = "false"


suffix = f"_hb" + args.output_suffix
command=f"{rosetta} " \
	f"-set_pdb_author 'Ryan Kibler (rdkibler@gmail.com)' " \
	f"-corrections::beta_nov16 true " \
	f"-holes:dalphaball /software/rosetta/DAlphaBall.gcc " \
	f"-out::file::pdb_comments " \
	f"-parser:protocol {protocol} " \
	f"-parser:script_vars {format_script_vars(script_vars_dict)} " \
	f"{input_lines} " \
	f"-nstruct 1 " \
	f"-overwrite 0 " \
	f"-out:file:scorefile_format json " \
	f"-no_nstruct_label " \
	f"-out:suffix {suffix} " \
	f"-out::path::all {outpath} " \
	f"{extra_lines} " \
	#f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"
	#f"-mute all " \

	#f"-unmute protocols.hbnet.HBNet " \


#actual run:
outfile = outpath + "/" + name + suffix + ".log"



run_command(command, outfile,to_file=False)
#run_command(command, outfile,to_file=args.debug)

# try:
# 	print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
# except KeyError:
# 	print("done %s" % (os.environ["SLURM_JOB_ID"]))




