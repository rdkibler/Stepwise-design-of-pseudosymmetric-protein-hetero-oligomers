#!/home/rdkibler/.conda/envs/pyro/bin/python3
import argparse

print("begin")

parser = argparse.ArgumentParser()
parser.add_argument("pdb")
parser.add_argument("--output_dir",type=str,help="default: output/",default="output/")

args = parser.parse_args()

print("args parsed")

import json

import glob
import subprocess
import pymol
pymol.finish_launching(['pymol','-qck']) #quiet, commandline, no plugins/pymolrc
import sys
sys.path.append("/home/rdkibler/scripts/pymol/will_pymol2")
import symgen
import tempfile
import os
from pathlib import Path
Path(args.output_dir).mkdir(parents=True, exist_ok=True)
name=args.pdb.split("/")[-1].split(".pdb")[0]

print("imports done")


pymol.cmd.load(args.pdb, "symmetric_design")
pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C4_A.pdb","original_chain_A")

pymol.cmd.copy_to("design_chain_A","symmetric_design and chain A")

pymol.cmd.select("sele1","design_chain_A and resi 201-276")
pymol.cmd.select("sele2","original_chain_A and resi 201-276")

pymol.cmd.align("polymer and name CA and sele2","polymer and name CA and sele1",quiet=0,object="aln",reset=1)

pymol.cmd.remove("design_chain_A and resi 240-276")
pymol.cmd.copy_to("design_chain_A","original and resi 240-314")
pymol.cmd.alter("design_chain_A","segi=''")
pymol.cmd.alter("design_chain_A","chain='A'")

#no need to renumber

#symgen.makecx(sel="design",n=4)


pymol.cmd.save("sesh.pse")


#pdbstr = pymol.cmd.get_pdbstr("SYM_design")
pdbstr = pymol.cmd.get_pdbstr("design_chain_A")
pymol.cmd.delete("*") #save some mem

print("assembly done")

import pyrosetta

pyrosetta.init("--in:file:fullatom -corrections::beta_nov16")

pose = pyrosetta.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose,pdbstr)

#pose.dump_pdb("asym.pdb.gz")

setup_for_symmetry = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover()
setup_for_symmetry.process_symmdef_file("/software/rosetta/latest/database/symmetry/cyclic/C4_Z.sym")
setup_for_symmetry.apply(pose)

pyrosetta.rosetta.core.pose.renumber_pdbinfo_based_on_conf_chains(pose)

print("loaded into pyrosetta")

#pose.dump_pdb("sym.pdb.gz")

extract_asymmetric_unit = pyrosetta.rosetta.protocols.symmetry.ExtractAsymmetricUnitMover()


xml_obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
<SCOREFXNS>
        <ScoreFunction name="hard_cart_cst" weights="beta_nov16_cart">
		<Reweight scoretype="coordinate_constraint" weight="1.0"/>
                <Reweight scoretype="cart_bonded" weight="1.5"/>
                <Reweight scoretype="pro_close" weight="0.0"/> #avoids double counting proline
		</ScoreFunction>
	<ScoreFunction name="hard" weights="beta_nov16"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
	<True name="all"/>
	<Index name="minimize_me_asym" resnums="204-276"/>
	<SymmetricalResidue name="minimize_me_sym" selector="minimize_me_asym"/>
	<Not name="lock_selector" selector="minimize_me_sym"/>
    </RESIDUE_SELECTORS>

    <TASKOPERATIONS>
	<InitializeFromCommandline name="init"/>
	<IncludeCurrent name="current"/>
	<LimitAromaChi2 name="arochi"/>
	<ExtraRotamersGeneric ex1="1" ex2="1" name="ex1_ex2"/>
	<OperateOnResidueSubset name="lock" selector="lock_selector">
		<PreventRepackingRLT/>
	</OperateOnResidueSubset>
    </TASKOPERATIONS>

    <MOVERS>
	<AddConstraints name="add_coord_csts">
        	<CoordinateConstraintGenerator
                	    name="local_movements_only"
	                    sd="0.6"
        	            residue_selector="all"/>
        </AddConstraints>

	    <ClearConstraintsMover name="remove_constraints" />

	#relax first with constraints on to work out any issues with the fusion we just made
    	SymMinMover name="cart_min" bondangle="1" bondlength="1" cartesian="1" scorefxn="hard_cart_cst" tolerance="0.0001" max_iter="2000" type="lbfgs_armijo_nonmonotone"  bb="1" chi="1" jump="0"/>

	#then

        <FastRelax name="relax" cartesian="1" scorefxn="hard_cart_cst" repeats="1" task_operations="init,current,arochi,ex1_ex2,lock">
		<MoveMap name="mm" bb="0" chi="0" jump="0">
			<ResidueSelector selector="minimize_me_sym" chi="1" bb="1" bondangle="1" bondlength="1"/>
		</MoveMap>
	</FastRelax>



    </MOVERS>
""")
add_coord_csts = xml_obj.get_mover("add_coord_csts")
relax = xml_obj.get_mover("relax")

print("parsed xml")

#Minimize the backbone with constraints to work out any issues from the fusion we just made
add_coord_csts.apply(pose)

relax.apply(pose)
print('relaxed')

pyrosetta.dump_pdb(pose, args.output_dir + "/" + name + "_rearmed.pdb.gz") 
