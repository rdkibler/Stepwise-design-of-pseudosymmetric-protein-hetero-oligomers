#!/usr/bin/python3

import uuid
import sys
import subprocess
import os
import re
import shlex

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

name = os.path.splitext(os.path.basename(input_pdb))[0]



outpath = os.path.abspath(f"output/{out_dir_name}")

from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)

protocol="xml/design_hbnets.xml"

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
with open(input_pdb,'r') as pdb:
	for line in pdb.readlines():
		if line[1:9] == "AtomPair":
			constraints.append(line[1:])
constraints_str = "".join(constraints)
#print(constraints_str)





#The CST files Scott writes are actually not functional (yay) so we gotta make a temp file
import tempfile
with tempfile.NamedTemporaryFile() as temp:
	temp.seek(0)
	temp.write(constraints_str.encode('utf-8'))
	temp.flush()

	native_pdb = "../../all_scaffolds/churro_9x25GB28GB_2_works.pdb"

	suffix = f"_designed"
	command=f"/home/bcov/ppi/builds/19-09-18/main/source/cmake/build_release_omp_hdf5/rosetta_scripts " \
		f"-corrections::beta_nov16 true " \
		f"-holes:dalphaball /home/bcov/ppi/tutorial_build/main/source/external/DAlpahBall/DAlphaBall.gcc " \
		f"-indexed_structure_store:fragment_store /home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h5 " \
		f"-out::file::pdb_comments " \
		f"-parser:protocol {protocol} " \
		f"-s {input_pdb} " \
		f"-in:file:native {native_pdb} " \
		f"-parser:script_vars cst_file={temp.name} " \
		f"-nstruct 1 " \
		f"-overwrite 1 " \
		f"-out:file:scorefile_format json " \
		f"-unmute all " \
		f"-mute all " \
		f"-no_nstruct_label " \
		f"-out:suffix {suffix} " \
		f"-out:file:scorefile {name}.sc " \
		f"-out::path::all {outpath} " \
		f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"



	#actual run:
	outfile = outpath + "/" + name + suffix + ".log"
	run_command(command, outfile)

	# try:
	# 	print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
	# except KeyError:
	# 	print("done %s" % (os.environ["SLURM_JOB_ID"]))




