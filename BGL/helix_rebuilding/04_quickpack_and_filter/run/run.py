#!/home/rdkibler/.conda/envs/pyro/bin/python3.7
import uuid
import sys
import subprocess
import os
import re
import shlex

input_pdb = sys.argv[1]
bin = sys.argv[2]

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

outpath = os.path.abspath(f"output/{bin}")

from pathlib import Path
Path(outpath).mkdir(parents=True, exist_ok=True)

protocol="xml/quikpak_n_filt.xml"

suffix = f"_hbfilt"

constraints=[]
with open(input_pdb,'r') as pdb:
        for line in pdb.readlines():
                if line[1:9] == "AtomPair":
                        constraints.append(line[1:])
constraints_str = "".join(constraints)

#The CST files Scott writes are actually not functional (yay) so we gotta make a temp file
import tempfile
with tempfile.NamedTemporaryFile(dir="/net/scratch/rdkibler/tmp") as temp:
	temp.seek(0)
	temp.write(constraints_str.encode('utf-8'))
	temp.flush()

	command=f"/home/bcov/ppi/builds/19-09-18/main/source/cmake/build_release_omp_hdf5/rosetta_scripts " \
		f"-corrections::beta_nov16 true " \
		f"-holes:dalphaball /home/bcov/ppi/tutorial_build/main/source/external/DAlpahBall/DAlphaBall.gcc " \
		f"-indexed_structure_store:fragment_store /home/brunette/DBs/hdf5/ss_grouped_vall_helix_shortLoop.h5 " \
		f"-parser:script_vars cst_file={temp.name} native={input_pdb}" \
		f"-out::file::pdb_comments " \
		f"-parser:protocol {protocol} " \
		f"-s {input_pdb} " \
		f"-in:file:native {input_pdb} " \
		f"-nstruct 1 " \
		f"-overwrite 0 " \
		f"-out:file:scorefile_format json " \
		f"-out:file:scorefile {name}.sc " \
		f"-mute all " \
		f"-no_nstruct_label " \
		f"-out:suffix {suffix} " \
		f"-out::path::all {outpath} " \
		f"-mute protocols.simple_moves.ScoreMover core.select.residue_selector.SecondaryStructureSelector core.select.residue_selector.PrimarySequenceNeighborhoodSelector"

	#actual run:
	outfile = outpath + "/" + name + suffix + ".log"
	run_command(command, outfile)

# try:
#       print("done %s" % (os.environ["SLURM_ARRAY_JOB_ID"] + "_" + os.environ["SLURM_ARRAY_TASK_ID"]))
# except KeyError:
#       print("done %s" % (os.environ["SLURM_JOB_ID"]))