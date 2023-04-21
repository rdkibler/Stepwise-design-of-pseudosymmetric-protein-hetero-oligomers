#!/home/rdkibler/.conda/envs/pyro/bin/python3

import argparse

parser = argparse.ArgumentParser("A script for performing recombination and resurfacing of any number of symmetric homologs")

parser.add_argument('homologs',nargs='+',help="A list of filenames in the order you want them to appear paired with residue numbers to use as cut points. Structures should have all symmetric chains. Example: RK388.pdb 122 RK390.pdb 122 RK410.pdb 132")
parser.add_argument('--recombine_only', action="store_true",default=False,help="only perform the recombination and skip the resurfacing. Useful for debugging")
parser.add_argument('--counter_clockwise', action="store_true",default=False,help="switch the order that chains are read in")
parser.add_argument('--output',type=str,default='output.pdb',help='filename of finished product')
parser.add_argument('--cyclic_symmetry',type=int,default=1,help="type of symmetry resulting fusion will have. REAL SYMMETRY, not pseudosymmetry. Use this only when identical chains are repeated (e.g. [AB]2)")
parser.add_argument('--debug_mode',action="store_true",default=False,help="prints a pymol session file to help show how the assembly occurred")

args = parser.parse_args()


if len(args.homologs) % 2 != 0:
	exit("Error! The list of homologs and cutpoints is not even numbered in length. Did you also include the cutpoints?")

filenames = []
interface_id_cutpoint_pairs = []
for i in range(int(len(args.homologs)/2)):
	filename_idx = i * 2
	cutpoint_idx = i * 2 + 1
	filename = args.homologs[filename_idx]
	try:
		cutpoint = int(args.homologs[cutpoint_idx])
	except:
		exit("Error! cutpoint \"{args.homologs[cutpoint_idx]}\" could not be cast to an int.")
	filenames.append(filename)
	interface_id_cutpoint_pairs.append([f"interface_{i}",cutpoint])
	#filename_cutpoint_pairs.append([filename,cutpoint])


import pymol
from pymol import cmd
import sys
import subprocess







def renumber(sele="all"):
	pdbstr = cmd.get_pdbstr(sele)
	process = subprocess.Popen(['perl',"/home/rdkibler/scripts/renum.pl"], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
	fixed_pdbstr = process.communicate(input=pdbstr.encode())[0].decode('utf-8')
	cmd.delete(sele)
	cmd.read_pdbstr(fixed_pdbstr,sele)




sys.path.insert(0,'/home/rdkibler/scripts/pymol/will_pymol2/')
from symgen import makecx
from pymol_util import alphabet
#/home/sheffler/pymol3/xyzMath.py
#run /home/sheffler/pymol3/xyzGeom.py
#run /home/sheffler/pymol3/pymol_util.py
#run /home/sheffler/pymol3/sym_util.py
#run /home/sheffler/pymol3/symgen.py

#makecx(sel='all', name=None, n=5)


num_unique_chains = len(filenames)
num_chains_list = []
#one by one, load them into pymol and do some preprocessing
for filename,(interface_id,cutpoint) in zip(filenames,interface_id_cutpoint_pairs):
	print(filename,interface_id,cutpoint)
	cmd.load(filename,interface_id+"_prime")
	num_chains_prime = len(cmd.get_chains(interface_id+"_prime")) #get_chains returns a list of chain letters
	num_chains_list.append(num_chains_prime)
	cmd.remove(f"{interface_id}_prime and not chain A")
	makecx(sel=interface_id+"_prime",name=interface_id,n=num_chains_prime) #just fix residue numbering so I can align on the cutpoints and only need to select chain numbers (b/c equivalent positions will have the same resnumber)


#check that all inputs are similar number of chains
num_chains_init = num_chains_list[0]
for num_chains_check in num_chains_list[1:]:
	if num_chains_check != num_chains_init:
		print("ERROR! not all inputs have the same number of chains. Exiting...")
		exit()


if args.counter_clockwise:
	"""
	if 4 chains...
	D stays the same
	C becomes A
	B becomes B (redundant but ok)
	X becomes C
	"""
	#print("A becomes X")
	cmd.alter('chain A',"chain='X'")
	for i in range(num_chains_init-1):
		cmd.alter(f"chain {alphabet[num_chains_init-i-1]}", f"chain='{alphabet[i]}'")
	#	print(alphabet[num_chains-i-1], "becomes", alphabet[i])
	#print("X becomes", alphabet[num_chains-1])
	cmd.alter("chain X",f"chain='{alphabet[num_chains_init-1]}'")

#perform alignment
print(f"interface_0 and resi {interface_id_cutpoint_pairs[0][1]}")
cmd.select("target",f"interface_0 and resi {interface_id_cutpoint_pairs[0][1]}")
for interface_id,cutpoint in interface_id_cutpoint_pairs[1:]:
	if args.debug_mode:
		cmd.create(f"{interface_id}_prealigned",interface_id)
	cmd.select("mobile",f"{interface_id} and resi {cutpoint}")
	cmd.align("mobile","target")

	#cmd.align(f"interface_{i+1} and resi {cutpoint}",f"interface_0 and resi {filename_cutpoint_pairs[0][1]}")
	#cmd.align(f"interface_{i+1}///{cutpoint}/CA",f"interface_0///{filename_cutpoint_pairs[0][1]}/CA")


#now extract fragments and form hybrid chain objects
chains = []
#want to copy N-terminus then C-terminus
#loop through chains in order to add the N-termini
for i,(interface_id,cutpoint) in enumerate(interface_id_cutpoint_pairs):
	cmd.create(f"chain_{alphabet[i]}",f"{interface_id} and resi 1-{cutpoint} and chain {alphabet[i]}")
	#renumber(f"chain_{alphabet[i]}")

	cmd.create(f"{interface_id}_Nfrag_{alphabet[i]}",f"{interface_id} and resi 1-{cutpoint} and chain {alphabet[i]}")
#now to C-terms
for i,(interface_id,cutpoint) in enumerate([interface_id_cutpoint_pairs[-1]] + interface_id_cutpoint_pairs[:-1]):
	#extract first, renumber the extraction, then join?
	cmd.create(f"{interface_id}_Cfrag_{alphabet[i]}",f"{interface_id} and resi {cutpoint+1}-100000 and chain {alphabet[i]}")
	#renumber(f"{interface_id}_Cfrag_{alphabet[i]}",start=cmd.count_atoms(f"chain_{alphabet[i]} and n. CA")+1,quiet=0)
	#cmd.copy_to(f"chain_{alphabet[i]}",f"{interface_id} and resi {cutpoint+1}-100000 and chain {alphabet[i]}")
	cmd.copy_to(f"chain_{alphabet[i]}",f"{interface_id}_Cfrag_{alphabet[i]}")
	renumber(f"chain_{alphabet[i]}") #yet this still doesn't seem to do the right thing??

	#correct chain id/segid so renumber works correctly. Shouldn't be necessary to do renumber here since I renumbered the fragment, but you know whatever. It needs to happen eventually for merging into 1 object
	chain_letter = alphabet[i]
	cmd.alter(f"chain_{chain_letter}",f"chain='{chain_letter}'")
	cmd.alter(f"chain_{chain_letter}",f"segi=''")




cmd.create("final",f"chain_A")
#cmd.delete(f"chain_A")
for chain_letter in alphabet[1:num_unique_chains]:
	cmd.copy_to("final",f"chain_{chain_letter}")
#	cmd.delete(f"chain_{chain_letter}")
cmd.alter("final","segi=''")
cmd.sort()

if args.debug_mode:
	cmd.save("assembly.pse")



if args.recombine_only:
	cmd.save(args.output,"final")
	exit("Exiting early because --recombine_only was set")



#now we go to pyrosetta
import pyrosetta
pyrosetta.init("-corrections::beta_nov16 true -out::file::pdb_comments ")
pose = pyrosetta.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, cmd.get_pdbstr("final"))

#cmd.delete("final")
cmd.delete("*")

if args.cyclic_symmetry > 1:
	print(f"applying C{args.cyclic_symmetry} symmetry")
	#will's makecx is only good for single chain jobs, so we gotta use rosetta
	#makecx(sel="final",n=args.cyclic_symmetry)
	sfsm = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover()
	sfsm.process_symmdef_file(f"/software/rosetta/latest/database/symmetry/cyclic/C{args.cyclic_symmetry}_Z.sym")
	sfsm.apply(pose)

chain_selectors = ""
mover_definitions = ""
mover_adds = ""

asym_chains = ",".join([f"ch{chain_letter}" for chain_letter in alphabet[:num_unique_chains]])
for chain in alphabet[:num_unique_chains]:
	chain_selectors += f"<Chain name=\"ch{chain}\" chains=\"{chain}\"/>"
	mover_definitions += f"<AddNetChargeConstraintMover name=\"net_charge_{chain}\" filename=\"/home/rdkibler/scripts/resurface/under_minus_two.charge\" selector=\"ch{chain}\" />"
	mover_adds += f"<Add mover=\"net_charge_{chain}\"/>\n"

import itertools
list_of_interfaces = []
interface_selectors = ""
for c1,c2 in itertools.combinations(alphabet[:num_unique_chains],2):
	interface_asym = f"interface{c1}{c2}_asym"
	interface = f"interface{c1}{c2}"
	list_of_interfaces.append(interface)
	interface_selectors += f"<InterfaceByVector name=\"{interface_asym}\" grp1_selector=\"ch{c1}\" grp2_selector=\"ch{c2}\" />"
	interface_selectors += f"<SymmetricalResidue name=\"{interface}\" selector=\"{interface_asym}\"/>"
list_of_interfaces_str = ",".join(list_of_interfaces)
join_interfaces = f"<Or name=\"interfaces\" selectors=\"{list_of_interfaces_str}\"/>"
#interface_selectors = f"<Or name='interface' selectors='interfaceAB,interfaceBC,interfaceCA'/>"

xmlstring = f"""
		<SCOREFXNS>

				<ScoreFunction name="sfxn" weights="beta_nov16_cart"/> # score function without weights for output

				<ScoreFunction name="up_ele" weights="beta_nov16"> # for designing surface
						<Reweight scoretype="fa_elec" weight="1.4"/>
						<Reweight scoretype="hbond_sc" weight="2.0" />
						<Reweight scoretype="res_type_constraint" weight="1.5"/> # For StructProfileMover
						<Reweight scoretype="netcharge" weight="1.0" />
						<Reweight scoretype="aa_composition" weight="1.0" /> # If you do use comp as well as netcharge, you'll need to change the comp file.
				</ScoreFunction>



		</SCOREFXNS>
		#############################


# Borrowed from ilutz [last revision 200407] on 200907
#modified 201120. Restrict packing around interface

		#############################
		# The score function section defines scorefunctions that will be used in Filters and Movers.
		# This can be used to define any of the scores defined in the path/to/rosetta/main/database.







		#############################
		# ResidueSelectors are used by movers, filters and task operations to dynamically select residues at run-time.
		# They are used to specify sets of residues based on multiple different properties.
		<RESIDUE_SELECTORS>

				# The LayerSelector lets a user select residues by burial.
				# Burial can be assessed by number of sidechain neighbors within a cone along the CA-CB vector (the default method), or by SASA


				Isaac's definitions
				<Layer name="surface" select_core="false" select_boundary="false" select_surface="true" use_sidechain_neighbors="true" core_cutoff="4.2" surface_cutoff="1.8"/>
				<Layer name="boundary" select_core="false" select_boundary="true" select_surface="false" use_sidechain_neighbors="true" core_cutoff="4.2" surface_cutoff="1.8"/>
				<Layer name="core" select_core="true" select_boundary="false" select_surface="false" use_sidechain_neighbors="true" core_cutoff="4.2" surface_cutoff="1.8"/>

				<Not name="not_surface" selector="surface" />
				<Not name="not_boundary" selector="boundary" />
				<Not name="not_core" selector="core" />


				# SecondaryStructureSelector selects all residues with given secondary structure
				<SecondaryStructure name="entire_helix" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="H" />
				<SecondaryStructure name="sheet" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="E" />
				<SecondaryStructure name="entire_loop" overlap="0" minH="3" minE="2" include_terminal_loops="true" use_dssp="true" ss="L" />

				<And name="helix_cap" selectors="entire_loop"> # Define a helix cap selection for layer design (C-term)
						<PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/>
				</And>
				<And name="helix_start" selectors="entire_helix"> # Define a helix start selection for layer design (N-term)
						<PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/>
				</And>
				<And name="helix" selectors="entire_helix"> # Define helix without its caps
						<Not selector="helix_start"/>
				</And>
				<And name="loop" selectors="entire_loop"> # Define loop without helix caps
						<Not selector="helix_cap"/>
				</And>

				<Neighborhood name="around_loop" distance="8.0" selector="entire_loop"/>
				<Not name="not_around_loop" selector="around_loop"/>

				<Or name="loop_and_around" selectors="loop,around_loop"/>
				<Not name="not_loop_and_around" selector="loop_and_around"/>

				<Or name="loop_and_caps" selectors="loop,helix_start,helix_cap"/>
		{chain_selectors}
		{interface_selectors}
		{join_interfaces}
		<Or name="asym_chains" selectors="{asym_chains}"/>
		<And name="asym_interfaces" selectors="asym_chains,interfaces"/>
		<Neighborhood name="generous_interfaces_asym" distance="10.0" selector="asym_interfaces"/>
		<SymmetricalResidue name="generous_interfaces" selector="generous_interfaces_asym"/>
				<ResiduePropertySelector name="polar" properties="POLAR"/>
				<Not name="non_polar" selector="polar"/>
				<Or name="pack_selector" selectors="boundary,generous_interfaces,non_polar"/>
				<Or name="lock_selector" selectors="core"/>


				<ResidueName name="pro_and_gly_positions" residue_name3="PRO,GLY" />

		</RESIDUE_SELECTORS>
		###############################


		<RESIDUE_LEVEL_TASK_OPERATIONS>
		PreventRepackingRLT name="PreventRepacking" />
		RestrictToRepackingRLT name="RestrictToRepacking" />

		# this mess is to add the current residue to the allowed aas list dynamically so that "no change" is always an option.
		# This should only be used on previously designed backbones (no poly X)
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT1"   aas="DEHKPQR"            />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT2"   aas="EHKQR"              />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT3"   aas="EHKNQRST"           />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT4"   aas="DEGHKNPQRST"        />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT5"   aas="ADEHIKLNPQRSTVWY"   />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT6"   aas="ADEHIKLNQRSTVWYM"   />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT7"   aas="DEFHIKLNQRSTVWY"    />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT8"   aas="ADEFGHIKLNPQRSTVWY" />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT9"   aas="AFILVWYP"           />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT10"  aas="AFILVWYM"      />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT11"  aas="FILVWY"        />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT12"  aas="AFGILPVWYM"    />
		<RestrictAbsentCanonicalAASExceptNativeRLT  name="RestrictAbsentCanonicalAASExceptNativeRLT13"  aas="DNSTP"              />
		</RESIDUE_LEVEL_TASK_OPERATIONS>


		###############################
		# TaskOperations are used by movers to tell the "packer" which residues/rotamers to use in reorganizing/mutating sidechains.
		# When used by certain Movers, the TaskOperations control what happens during packing, usually by restriction "masks".
		# TaskOperations can also be used by movers to specify sets of residues to act upon in non-packer contexts.
		<TASKOPERATIONS>
				<IncludeCurrent name="current"/> # Tell the packer to also consider the input rotamer.
				<LimitAromaChi2 name="arochi" /> # Prevents use the rotamers of PHE, TYR and HIS that have chi2 far from 90. These rotamers are acceptable to the Dunbrack energy term due to smoothing of the energy landscape but not actually physically common without some other rearrangements Rosetta fails to capture.
				<ExtraRotamersGeneric name="ex1_ex2" ex1="1" ex2="1"/> # Add increased sampling for chi1 and chi2 rotamers

				<DesignRestrictions name="layer_design">
					<Action selector_logic="surface AND helix_start"    aas="DEHKPQR"               residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT1"/>
					<Action selector_logic="surface AND helix"          aas="EHKQR"                 residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT2"/>
					<Action selector_logic="surface AND sheet"          aas="EHKNQRST"              residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT3"/>
					<Action selector_logic="surface AND loop"           aas="DEGHKNPQRST"           residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT4"/>
					<Action selector_logic="boundary AND helix_start"   aas="ADEHIKLNPQRSTVWY"      residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT5"/>
					<Action selector_logic="boundary AND helix"         aas="ADEHIKLNQRSTVWYM"      residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT6"/>
					<Action selector_logic="boundary AND sheet"         aas="DEFHIKLNQRSTVWY"       residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT7"/>
					<Action selector_logic="boundary AND loop"          aas="ADEFGHIKLNPQRSTVWY"    residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT8"/>
					<Action selector_logic="core AND helix_start"       aas="AFILVWYP"              residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT9"/>
					<Action selector_logic="core AND helix"             aas="AFILVWYM"              residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT10"/>
					<Action selector_logic="core AND sheet"             aas="FILVWY"                residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT11"/>
					<Action selector_logic="core AND loop"              aas="AFGILPVWYMN"           residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT12"/>
					<Action selector_logic="helix_cap"                  aas="DNSTP"                 residue_level_operations="RestrictAbsentCanonicalAASExceptNativeRLT13"/>
				</DesignRestrictions>

				<ConsensusLoopDesign name="disallow_non_abego_aas"/> # For loop design

				<OperateOnResidueSubset name="pack" selector="pack_selector">
			<RestrictToRepackingRLT/>
		</OperateOnResidueSubset>

				<OperateOnResidueSubset name="lock" selector="lock_selector">
			<PreventRepackingRLT/>
		</OperateOnResidueSubset>

				<OperateOnResidueSubset name="restrict_PRO_GLY" selector="pro_and_gly_positions">
						<RestrictToRepackingRLT/>
				</OperateOnResidueSubset>


		</TASKOPERATIONS>


		###############################
		# Define the stuff that will do stuff to the protein (e.g. changes sequence, geometry, etc...).
		<MOVERS>

				### repeat this for as many chains as you have
		{mover_definitions}
				<PackRotamersMover name="redesign_surface" scorefxn="up_ele" task_operations="current,arochi,layer_design,ex1_ex2,restrict_PRO_GLY,pack,lock,disallow_non_abego_aas" />

				FastDesign name="redesign_surface" scorefxn="up_ele" task_operations="current,arochi,layer_design,ex1_ex2,restrict_PRO_GLY,pack,lock,disallow_non_abego_aas" repeats="1" relaxscript="MonomerDesign2019" >
					MoveMap name="what_moves" bb="false" chi="true" jump="false" />
				/FastDesign>

				<MinMover name="min_all" scorefxn="sfxn" chi="true" bb="true" jump="ALL" cartesian="true" type="lbfgs_armijo_nonmonotone" tolerance="0.0001" max_iter="200" />

		</MOVERS>


	<PROTOCOLS>
		#cart min to fix bb after fusion
				<Add mover="min_all"/>
		#add a AddNetChargeConstraintMover for each chain in pose
		{mover_adds}
				<Add mover="redesign_surface"/>
		#another cart min to rescore and make everything nice and perfect before outputting
				<Add mover="min_all"/>
		#min, repack, min is like poor man's fast design

	</PROTOCOLS>


"""
print(xmlstring)
xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xmlstring)


protocol = xml.get_mover("ParsedProtocol")
protocol.apply(pose)
pose.dump_pdb(args.output)

