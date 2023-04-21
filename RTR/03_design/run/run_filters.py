#!/usr/bin/python3

import uuid
import sys
import subprocess
import os
import re
import shlex
import gzip
import json

input_pdb = sys.argv[1]
out_dir_name = sys.argv[2]


production = False

def run_command(command, outfile):
	if production:
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
			result = subprocess.run(shlex.split(command), stdout=out, stderr=out)
			#result = subprocess.run(shlex.split(command), stdout=subprocess.PIPE)
			#out.write(result.stdout.decode('utf-8'))
	else:
		subprocess.run(shlex.split(command))

#all the outputs are "pdb"
name = os.path.basename(input_pdb).split(".pdb")[0]


outpath = os.path.abspath(f"/net/scratch/rdkibler/210126_nmp_hbnet_packing/{out_dir_name}")

from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)

protocol="xml/filtering_only.xml"

import glob
import os

net_num = input_pdb.split("_")[-1].split(".")[0].lstrip("0")
#try :
#	cst_file = glob.glob(f"{os.path.dirname(input_pdb)}/*_hb_network_{net_num}.cst")[0]
#except:
#	cst_file = f"{os.path.splitext(input_pdb)[0]}_hb_network_1.cst"


#these constraints are written to the pdb so I can actually turn off the cst write!
#now I can pre-filter
constraints=[]
HBNet_resis = []
with open(input_pdb,'r') as pdb:
	for line in pdb.readlines():
		if line[1:9] == "AtomPair":
			constraints.append(line[1:])
		if " HBNet" in line:
			HBNet_resis.append(line.split()[2])
constraints_str = "".join(constraints)
#print(constraints_str)

#descor = "/home/rdkibler/projects/tiara_gen2/toroids/05_2_ctor-24x_parametric/03_generate_hbnets/design_correspondance.json"
#with open(descor,'r') as f:
#	connection_lookup = json.load(f)

#if connection_lookup[os.path.basename(input_pdb)].split("/")[6].split("_")[-1] == "C":
#	terminal_helix_selection_start='-2,H,M'
#	terminal_helix_selection_stop='c_term'
#elif connection_lookup[os.path.basename(input_pdb)].split("/")[6].split("_")[-1] == "N":
#	terminal_helix_selection_start='n_term'
#	terminal_helix_selection_stop='2,H,M'
#else:
#	print("NOT RECOGNIZED:")
#	print(connection_lookup[os.path.basename(input_pdb)].split("/")[6].split("_")[-1])

#print("might need reworking")
#exit()


#The CST files Scott writes are actually not functional (yay) so we gotta make a temp file
import tempfile
with tempfile.NamedTemporaryFile() as temp:
	temp.seek(0)
	temp.write(constraints_str.encode('utf-8'))
	temp.flush()

	suffix = f"_designed"
	#command=f"/home/bcov/ppi/builds/19-09-18/main/source/cmake/build_release_omp_hdf5/rosetta_scripts " \
	#	f"-holes:dalphaball /home/bcov/ppi/tutorial_build/main/source/external/DAlpahBall/DAlphaBall.gcc " \

#	command=f"/net/software/rosetta/versions/v2020.25-dev61318/main/source/build/src/release/linux/5.4/64/x86/gcc/9/hdf5/rosetta_scripts.hdf5.linuxgccrelease " \
	command=f"/net/software/rosetta/versions/v2020.50-dev61505/main/source/build/src/release/linux/5.4/64/x86/gcc/9/hdf5/rosetta_scripts.hdf5.linuxgccrelease " \
		f"-corrections::beta_nov16 true " \
		f"-indexed_structure_store:fragment_store /home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h5 " \
		f"-holes:dalphaball /software/rosetta/DAlphaBall.gcc " \
		f"-out::file::pdb_comments " \
		f"-parser:protocol {protocol} " \
		f"-s {input_pdb} " \
		f"-parser:script_vars cst_file={temp.name} HBNet_resnums={','.join(HBNet_resis)} " \
		f"-nstruct 1 " \
		f"-overwrite 0 " \
		f"-out:file:scorefile_format json " \
		f"-unmute all " \
		f"-no_nstruct_label " \
		f"-out:suffix {suffix} " \
		f"-out:file:scorefile {name}.sc " \
		f"-out::path::all {outpath} " \
		f"-out:pdb_gz 1 " \
		f"-mute all " \

		#f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"



	#actual run:
	outfile = outpath + "/" + name + suffix + ".log"
	run_command(command, outfile)

	# try:
	# 	print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
	# except KeyError:
	# 	print("done %s" % (os.environ["SLURM_JOB_ID"]))




