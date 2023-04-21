#!/usr/bin/python3

import uuid
import sys
import subprocess
import os
import re
import shlex
import pymol

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
outpath = os.path.abspath(f"/net/scratch/rdkibler/201223_rotor_nmp/n{n:05}_pert{pertscale}_{nmodes}modes/")
from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)

protocol="xml/perterb_interf.xml"
ft_file = input_pdb.replace(".pdb",".ft")

suffix = f"_nm_{n:05}_{pertscale}_{nmodes}"

pymol.cmd.load(input_pdb,'obj1')
seq = "".join([line.strip() for line in pymol.cmd.get_fastastr('obj1').split("\n") if ">" not in line])

print(seq)

#command=f"/home/bcov/ppi/builds/19-09-18/main/source/cmake/build_release_omp_hdf5/rosetta_scripts " \
#command=f"/software/rosetta/latest/bin/rosetta_scripts.hdf5.linuxgccdebug " \
#command=f"/mnt/home/rdkibler/projects/normal_modes_MPM/Rosetta/main/source/build/src/debug/linux/5.4/64/x86/gcc/7/hdf5/rosetta_scripts.hdf5.linuxgccdebug " \
command=f"/home/rdkibler/projects/normal_modes_MPM/Rosetta/main/source/bin/rosetta_scripts.hdf5.linuxgccrelease " \
	f"-corrections::beta_nov16 true " \
	f"-out::file::pdb_comments " \
	f"-parser:protocol {protocol} " \
	f"-parser:script_vars ft_file={ft_file} pertscale={pertscale} nmodes={nmodes} mm={mm} fasta={seq} " \
	f"-s {input_pdb} " \
	f"-in:file:native {input_pdb} " \
	f"-nstruct 1 " \
	f"-overwrite 0 " \
	f"-out:file:scorefile_format json " \
	f"-out:suffix {suffix} " \
	f"-out::path::all {outpath} " \

	#f"-mute all " \

	#f"-packing::unboundrot {input_pdb} " \

	#f"-out::level 9001 " \
	#f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"
	#f"-no_nstruct_label " \
	#f"-no_scorefile 1 " \




#actual run:
outfile = outpath + "/" + name + suffix + ".log"
run_command(command, outfile)

# try:
# 	print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
# except KeyError:
# 	print("done %s" % (os.environ["SLURM_JOB_ID"]))

