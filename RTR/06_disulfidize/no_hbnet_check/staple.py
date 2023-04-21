#!/home/rdkibler/.conda/envs/pyro/bin/python

import glob
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_dir")
parser.add_argument("output_dir")
args = parser.parse_args()

import json
import pandas as pd
import pyrosetta
pyrosetta.init()

from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector

from pyrosetta.rosetta.protocols.symmetry import DetectSymmetry
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

from stapler import *

# Preset ResidueSelectors
default_residue_selectors = [TrueResidueSelector(), TrueResidueSelector()]
interface_residue_selectors = [ChainSelector('A'), ChainSelector('B')]
interface_or_internal_residue_selectors = [ChainSelector('A'), ChainSelector('A,B')]
only_binder_residue_selectors = [ChainSelector('B'), ChainSelector('B')]
not_on_loops = [SecondaryStructureSelector('HE'), SecondaryStructureSelector('HE')]
not_on_loops_across_interface = [AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('A')),
                                 AndResidueSelector(SecondaryStructureSelector('HE'),ChainSelector('B'))]












for pdb in glob.glob(os.path.join(args.input_dir,'*.pdb*')):

	pdb_filename = os.path.basename(pdb).split(".pdb")[0]

	non_hbnet_chain_A = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector()
	non_hbnet_chain_B = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector()
	non_hbnet_chain_A.add_residue_selector(ChainSelector('A'))
	non_hbnet_chain_B.add_residue_selector(ChainSelector('B'))

	interface_residue_selectors_non_hbnet = [non_hbnet_chain_A,non_hbnet_chain_B]

	# Initialize the native disulfide stapler with defaults.
	native_disulfide_stapler = NativeDisulfideStapler(
	    residue_selectors=interface_residue_selectors_non_hbnet,
	    minimum_sequence_distance=0
	)


	pose = pyrosetta.pose_from_file(pdb)

	#Preset Movers for Symmetry
	DetectSymmetry().apply(pose)
	#SetupForSymmetryMover('C2.sym').apply(pose)

	for i, crosslinked_pose in enumerate(native_disulfide_stapler.apply(pose)):
		crosslinked_pose.dump_pdb(os.path.join(args.output_dir,f"{pdb_filename}_staple_{i}.pdb.gz"))
