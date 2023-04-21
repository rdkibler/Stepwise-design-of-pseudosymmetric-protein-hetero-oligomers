#!/home/rdkibler/.conda/envs/pyro/bin/python3
import pymol
import sys


import tempfile
import pyrosetta as pyro
import os
pyro.init("--in:file:fullatom --beta -corrections::beta_nov16")
'''
http://pymolwiki.org/index.php/Renumber

(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

def renumber(selection='all', start=1, startsele=None, quiet=1):
    '''
DESCRIPTION

    Set residue numbering (resi) based on connectivity.

ARGUMENTS

    selection = string: atom selection to renumber {default: all}

    start = integer: counting start {default: 1}

    startsele = string: residue to start counting from {default: first in
    selection}
    '''
    start, quiet = int(start), int(quiet)
    model = pymol.cmd.get_model(selection)
    pymol.cmd.iterate(selection, 'atom_it.__next__().model = model',
                space={'atom_it': iter(model.atom)})
    if startsele is not None:
        startidx = cmd.index('first (' + startsele + ')')[0]
        for atom in model.atom:
            if (atom.model, atom.index) == startidx:
                startatom = atom
                break
        else:
            print(' Error: startsele not in selection')
            raise pymol.CmdException
    else:
        startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [('SG', 'SG')]:
                traverse(other, resi)
    traverse(startatom, start)
    pymol.cmd.alter(selection, 'resi = atom_it.__next__().resi',
              space={'atom_it': iter(model.atom)})
    if not quiet:
        print(' Renumber: range (%d to %d)' % tuple(minmax))











input_obj = pymol.cmd.load(sys.argv[1], "input")

chain_A_start = 115 #chain A motif will alwasy be correct
chain_A_end = 171

prot_len = pymol.cmd.count_atoms("input and n. CA")
chA_len = pymol.cmd.count_atoms("input and chain A and n. CA")
chB_len = pymol.cmd.count_atoms("input and chain B and n. CA")

chain_B_start = prot_len - 227
chain_B_end = prot_len - 227 + 56

chain_A_obj = pymol.cmd.create("chA","input and chain A")
chain_B_obj = pymol.cmd.create("chB","input and chain B")

#pymol.cmd.save("testA1.pdb", "chA")
#pymol.cmd.save("testB1.pdb", "chB")

pymol.cmd.select("alnA", f"chA and resi {chain_A_start}-{chain_A_end}")
pymol.cmd.select("alnB", f"chB and resi {chain_B_start}-{chain_B_end}")

pymol.cmd.align("polymer and name CA and (alnA)","polymer and name CA and (alnB)",quiet=0,object="aln_alnA_to_alnB",reset=1)

#pymol.cmd.save("testA2.pdb", "chA")
#pymol.cmd.save("testB2.pdb", "chB")
#remove chB and resi 340-510
#remove chA and resi 1-114

#remove chA and resi 1-140 #half way
#remove chB and resi 309-510
pymol.cmd.remove(f"chA and resi 1-140") #this is ok
pymol.cmd.remove(f"chB and resi {prot_len - 201}-{prot_len}")

#pymol.cmd.save("testA3.pdb", "chA")
#pymol.cmd.save("testB3.pdb", "chB")

#pymol.cmd.create("merge","chB",zoom=0)
#pymol.cmd.copy_to("merge","chA",zoom=0,quiet=0)

pymol.cmd.select('C',f"n. c and chB and i. {prot_len - 202}",enable=1)
pymol.cmd.select('N',"chA`1",enable=1) #this always ok

#pymol.cmd.fuse("C","N",mode=0)
pymol.cmd.copy_to("merge","chB")
pymol.cmd.copy_to("merge","chA")
pymol.cmd.alter("merge","chain='A'")

pymol.cmd.sort()

renumber('merge')

#alter chA, chain="A"
pymol.cmd.alter("merge","segi=''")


#pymol.cmd.load("/home/rdkibler/projects/tiara_gen2/toroids/all_scaffolds/trimer/churro-9x_sub3.pdb","ref")
#pymol.cmd.align("polymer and name CA and (chA)","polymer and name CA and (ref and chain A)")
#pymol.cmd.save("temp.pdb", "merge")


file = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
file.close()

pymol.cmd.save(file.name, "merge")
pymol.cmd.delete("*") #save some mem

pose = pyro.pose_from_file(file.name)
setup_for_symmetry = pyro.rosetta.protocols.symmetry.SetupForSymmetryMover()
setup_for_symmetry.process_symmdef_file("C3_Z.sym")
setup_for_symmetry.apply(pose)

print("Don't do cartesian relax for sidechains. Only do it on backbone to fix it up. Then do normal relax on the sidechains and backbone")

xml_obj = pyro.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
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
    
    <MOVERS>
	    <AddConstraintsToCurrentConformationMover CA_only="true" coord_dev="2.0" name="add_coord_cst" use_distance_cst="0"/>
	    <ClearConstraintsMover name="remove_constraints" />

    	<SymMinMover name="cart_min" bondangle="1" bondlength="1" cartesian="1" scorefxn="hard_cart_cst" tolerance="0.0001" max_iter="2000" type="lbfgs_armijo_nonmonotone"  bb="1" chi="0" jump="ALL"/>

    	<SwitchResidueTypeSetMover name="fa" set="fa_standard"/>
        <FastRelax name="relax" scorefxn="hard" cartesian="false" repeats="1"/>
    </MOVERS>
 
""")
fa = xml_obj.get_mover("fa")
add_coord_cst = xml_obj.get_mover("add_coord_cst")
cart_min = xml_obj.get_mover("cart_min")
remove_constraints = xml_obj.get_mover("remove_constraints")
relax = xml_obj.get_mover("relax")

fa.apply(pose)
add_coord_cst.apply(pose)
cart_min.apply(pose)
remove_constraints.apply(pose)
relax.apply(pose)

pose.dump_pdb(os.path.splitext(os.path.basename(sys.argv[1]))[0] + "_sym.pdb")






os.unlink(file.name)