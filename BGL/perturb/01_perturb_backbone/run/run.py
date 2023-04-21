#!/usr/bin/python3

import uuid
import sys
import subprocess
import os
import re
import shlex

input_pdb = sys.argv[1]
n=int(sys.argv[2])
pertscale = sys.argv[3]
nmodes = sys.argv[4]
mm = sys.argv[5]

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

#outpath = os.path.abspath(f"output/{name}_{pertscale}_{nmodes}_{mm}")
#from pathlib import Path
#Path(outpath).mkdir(parents=True, exist_ok=True)
outpath = os.path.abspath("output/")
from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)

protocol="xml/perterb_interf.xml"
ft_file = input_pdb.replace(".pdb",".ft")

suffix = f"_nm"

command=f"/home/bcov/ppi/builds/19-09-18/main/source/cmake/build_release_omp_hdf5/rosetta_scripts " \
	f"-corrections::beta_nov16 true " \
	f"-holes:dalphaball /home/bcov/ppi/tutorial_build/main/source/external/DAlpahBall/DAlphaBall.gcc " \
	f"-indexed_structure_store:fragment_store /home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h5 " \
	f"-out::file::pdb_comments " \
	f"-parser:protocol {protocol} " \
	f"-parser:script_vars ft_file={ft_file} pertscale={pertscale} nmodes={nmodes} mm={mm} " \
	f"-s {input_pdb} " \
	f"-in:file:native {input_pdb} " \
	f"-nstruct 1 " \
	f"-no_nstruct_label " \
	f"-overwrite 0 " \
	f"-out:file:scorefile_format json " \
	f"-no_scorefile 1 " \
	f"-mute all " \
	f"-out:suffix {suffix}_{n:05} " \
	f"-out::path::all {outpath} " \
	f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"




#actual run:
outfile = outpath + "/" + name + suffix + ".log"
run_command(command, outfile)

# try:
# 	print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
# except KeyError:
# 	print("done %s" % (os.environ["SLURM_JOB_ID"]))




