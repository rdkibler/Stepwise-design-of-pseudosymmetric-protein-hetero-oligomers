#!/home/rdkibler/.conda/envs/pyro/bin/python3
import argparse

print("begin")

parser = argparse.ArgumentParser()
parser.add_argument("pdb")
parser.add_argument("Cx",type=int,choices=[1,2,3,4,6,8,12],help="Number of submits in final oligomer after inserting interface")
parser.add_argument("--hbnet_mapping",type=str,help="a json file which maps designed fiels to the input hbnet files for the sake a scraping cst info")
parser.add_argument("--output_dir",type=str,help="default: output/",default="output/")

args = parser.parse_args()

print("args parsed")

import json

import glob
import subprocess
import pymol
pymol.finish_launching(['pymol','-qck']) #quiet, commandline, no plugins/pymolrc
import sys
import tempfile
import pyrosetta as pyro
import os
from pathlib import Path
Path(args.output_dir).mkdir(parents=True, exist_ok=True)
name=args.pdb.split("/")[-1].split(".pdb")[0]

print("imports done")

if args.hbnet_mapping == None:
	print("performing mapping manually. It's not very efficient...")
	hbnet_paths = glob.glob("/net/scratch/rdkibler/201224_rotor_hbnet/*/*/*.pdb")
	hbnet_mapping = {path.split("/")[-1].split(".pdb")[0] + "_designed":path for path in hbnet_paths}
else:
	with open(args.hbnet_mapping) as f:
		hbnet_mapping = json.load(f)

print("mapping loaded")

hbnet_path = hbnet_mapping[name]

constraints = []
hbnet_lines = []
with open(hbnet_path, 'r') as hbp:
	for line in hbp:
		if line[1:9] == "AtomPair":
			constraints.append(line[1:])
		elif " HBNet" in line:
			hbnet_lines.append(line)

len_chain_A = 66
interface_N_homology_amt = 32
interface_N_homology_start = len_chain_A - interface_N_homology_amt
interface_N_homology_end = len_chain_A

interface_C_homology_start = len_chain_A + 1
interface_C_homology_amt = 29
interface_C_homology_end = interface_C_homology_start + interface_C_homology_amt

#convert constraints
adjusted_constraints = []
for constraint in constraints:
	constraint_parts = constraint.split()
	if int(constraint_parts[2]) <= len_chain_A:
		constraint_parts[2] = f"{constraint_parts[2]}A"
	else:
		constraint_parts[2] = f"{int(constraint_parts[2]) + len_chain_A}B"

	if int(constraint_parts[4]) <= len_chain_A:
		constraint_parts[4] = f"{constraint_parts[4]}A"
	else:
		constraint_parts[4] = f"{int(constraint_parts[4]) + len_chain_A}B"
	adjusted_constraints.append(" ".join(constraint_parts) + "\n")

#print("old")
#print("".join(constraints))
#print("new")
#print("".join(adjusted_constraints))


input_obj = pymol.cmd.load(args.pdb, "interface")


if args.Cx == 1:
	raise NotImplementedError()
elif args.Cx == 2:
	backbone_obj = pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C2_A.pdb", "bb")
elif args.Cx == 3:
	backbone_obj = pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C3_A.pdb", "bb")
elif args.Cx == 4:
	backbone_obj = pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C4_A.pdb", "bb")
elif args.Cx == 6:
	backbone_obj = pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C6_A.pdb", "bb")
elif args.Cx == 8:
	backbone_obj = pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/rotors/01_RC4_20_higher_order/01_single_chain_inputs/RC4_20_C8_A.pdb", "bb")
elif args.Cx == 12:
	raise NotImplementedError("this will take more work b/c this has no homology")


pymol.cmd.copy_to("interface_N","interface and chain A")
pymol.cmd.copy_to("interface_C","interface and chain B")


bb_N_homology_start = 34
bb_N_homology_end = 66

bb_len = pymol.cmd.count_atoms("bb and n. CA")
bb_C_homology_start = bb_len - 181
bb_C_homology_end = bb_len - 152

pymol.cmd.select("bb_N_homology",f"bb and i. {bb_N_homology_start}-{bb_N_homology_end}")
pymol.cmd.select("interface_N_homology",f"interface_N and i. {interface_N_homology_start}-{interface_N_homology_end}")
pymol.cmd.select("bb_C_homology",f"bb and i. {bb_C_homology_start}-{bb_C_homology_end}")
pymol.cmd.select("interface_C_homology",f"interface_C and i. {interface_C_homology_start}-{interface_C_homology_end}")

pymol.cmd.align("polymer and name CA and (interface_N_homology)","polymer and name CA and (bb_N_homology)",quiet=1,object="aln_N",reset=1)
pymol.cmd.align("polymer and name CA and (interface_C_homology)","polymer and name CA and (bb_C_homology)",quiet=1,object="aln_C",reset=1)

pymol.cmd.copy_to("assembled_asu","interface_N")
pymol.cmd.copy_to("assembled_asu",f"bb and i. {bb_N_homology_end + 1}-{bb_C_homology_start-1}")
pymol.cmd.copy_to("assembled_asu",f"interface_C")
pymol.cmd.alter("assembled_asu","chain='A'")
pymol.cmd.alter("assembled_asu","segi=''")





#file = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
#file.close()

choppy_pdbstr = pymol.cmd.get_pdbstr("assembled_asu")
pymol.cmd.delete("*") #save some mem

process = subprocess.Popen(['perl',"/home/rdkibler/scripts/renum.pl"], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
fixed_pdbstr = process.communicate(input=choppy_pdbstr.encode())[0].decode('utf-8')

print("assembly done")

pyro.init("--in:file:fullatom -corrections::beta_nov16 --mute all")
pose = pyro.Pose()
pyro.rosetta.core.import_pose.pose_from_pdbstring(pose,fixed_pdbstr)

setup_for_symmetry = pyro.rosetta.protocols.symmetry.SetupForSymmetryMover()
setup_for_symmetry.process_symmdef_file("/software/rosetta/latest/database/symmetry/cyclic/C4_Z.sym")
setup_for_symmetry.apply(pose)

pyro.rosetta.core.pose.renumber_pdbinfo_based_on_conf_chains(pose)

extract_asymmetric_unit = pyro.rosetta.protocols.symmetry.ExtractAsymmetricUnitMover()

cst_file = tempfile.NamedTemporaryFile(delete=False)
constraints_str = "".join(adjusted_constraints)
with open(cst_file.name,'w') as f:
	f.seek(0)
	f.write(constraints_str)
	f.flush()

pyro_csts = []
#print(pose)




#for line in adjusted_constraints:
#	print(line)
#	cst = pyro.rosetta.core.scoring.constraints.AtomPairConstraint()
#	cst.read_data(pyro.rosetta.std.istringstream(line))
#	print(cst.dist(pose))


print("loaded into pymol")

xml_obj = pyro.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(f"""
<SCOREFXNS>
        <ScoreFunction name="hard_cart_cst" weights="beta_nov16_cart">
		<Reweight scoretype="aa_composition" weight="1.0"/>
		<Reweight scoretype="coordinate_constraint" weight="1.0"/>
                <Reweight scoretype="atom_pair_constraint" weight="1.0"/>
                <Reweight scoretype="angle_constraint" weight="1.0"/>
                <Reweight scoretype="dihedral_constraint" weight="1.0"/>
                <Reweight scoretype="cart_bonded" weight="1.5"/>
                <Reweight scoretype="pro_close" weight="0.0"/> #avoids double counting proline
		</ScoreFunction>
		<ScoreFunction name="hard" weights="beta_nov16"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
        <True name="all"/>
    </RESIDUE_SELECTORS>

    <TASKOPERATIONS>
        <InitializeFromCommandline name="init"/>
        <IncludeCurrent name="current"/>
        <LimitAromaChi2 name="arochi"/>
        <ExtraRotamersGeneric ex1="1" ex2="1" name="ex1_ex2"/>
    </TASKOPERATIONS>

    <MOVERS>
	    <AddConstraintsToCurrentConformationMover CA_only="true" coord_dev="1.0" name="add_coord_cst" use_distance_cst="0"/>

		<AddConstraints name="add_coord_csts">
                <CoordinateConstraintGenerator
                    name="local_movements_only"
                    sd="0.6"
                    residue_selector="all"/>
        </AddConstraints>
        
        <AddConstraints name="add_hbnet_csts">
            <FileConstraintGenerator name="hbnet_csts" filename="{cst_file.name}"/>
        </AddConstraints>

	    <ClearConstraintsMover name="remove_constraints" />

        #relax first with constraints on to work out any issues with the fusion we just made
    	<SymMinMover name="cart_min" bondangle="1" bondlength="1" cartesian="1" scorefxn="hard_cart_cst" tolerance="0.0001" max_iter="2000" type="lbfgs_armijo_nonmonotone"  bb="1" chi="1" jump="0"/>
        #then
    	<SwitchResidueTypeSetMover name="fa" set="fa_standard"/>
        <FastRelax name="relax" scorefxn="hard" cartesian="false" repeats="1" task_operations="init,current,arochi,ex1_ex2"/>

    </MOVERS>
""")
add_coord_csts = xml_obj.get_mover("add_coord_csts")
add_hbnet_csts = xml_obj.get_mover("add_hbnet_csts")
cart_min = xml_obj.get_mover("cart_min")
remove_csts = xml_obj.get_mover("remove_constraints")
relax = xml_obj.get_mover("relax")

print("parsed xml")

#Minimize the backbone with constraints to work out any issues from the fusion we just made
add_hbnet_csts.apply(pose)

#Read constraints before minimization
pyro_csts = pose.constraint_set().get_all_constraints()
cst_data = {}
for cst in pyro_csts:
	cst_data[f"{cst.atom1().rsd()}_{cst.atom1().atomno()},{cst.atom2().rsd()}_{cst.atom2().atomno()},before"] = {"dist":cst.dist(pose), "score":cst.score(pose)}

add_coord_csts.apply(pose)
cart_min.apply(pose) #fix backbone bonds (hopefully)

#relax without constraints to check for packing of hbonds. If they move I guess the hbnets are underpacked
remove_csts.apply(pose) 
relax.apply(pose)

#now read out those csts
for cst in pyro_csts:
	cst_data[f"{cst.atom1().rsd()}_{cst.atom1().atomno()},{cst.atom2().rsd()}_{cst.atom2().atomno()},after"] = {"dist":cst.dist(pose), "score":cst.score(pose)}

pyro.dump_pdb(pose, args.output_dir + "/" + name + "_relax_sym.pdb") #DAN doesn't take pdb.gzs, so I'll just write these directly as pdbs

os.unlink(cst_file.name)

#write and do analysis later
with open(args.output_dir + "/" + name + "_cst_data.json", 'w') as f:
	f.write(json.dumps(cst_data)+"\n")

print("done")
