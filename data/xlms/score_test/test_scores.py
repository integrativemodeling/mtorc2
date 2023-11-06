from pathlib import Path
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.algebra
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
import sys
import math
import pandas as pd
sys.path.append("/home/matthew/mtorc2/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/src")
from restraints import get_xl_restraint


def create_rigid_body(
        mol,
        component,
        dof
):
    rbs = list()
    max_trans, max_rot = 0,0
    print(component, max_trans, max_rot)

    if component == "KINASE":
        rigid_parts_list = list()
        # The first entry is the N lobe and the second entry is the C lobe.
        for start, end in [("1", "89"), ("90", "335")]:
            rigid_parts_list.append(
                mol.residue_range(
                    a=start,
                    b=end
                )
            )

        non_rigid_parts_list = list()
        # The first entry is the short loop connecting the N and C lobes. The second entry is the tail region of the C lobe (which may or not be ommitted from the structure).
        for start, end in [("84", "89"), ("273", "335")]:
            non_rigid_parts_list.append(
                mol.residue_range(
                    a=start,
                    b=end
                )
            )
        # non_rigid_parts_list.append(None)

    else:
        rigid_parts_list = [mol.get_residues()]
        non_rigid_parts_list = [mol.get_non_atomic_residues()]

    for i in range(len(rigid_parts_list)):
        rigid_parts = rigid_parts_list[i]
        non_rigid_parts = non_rigid_parts_list[i]

        print("RIGID BODIES")
        print(len(rigid_parts))
        print(len(non_rigid_parts))

        rb_movers, rb = dof.create_rigid_body(
            rigid_parts=rigid_parts,
            nonrigid_parts=non_rigid_parts,
            max_trans=max_trans,
            max_rot=max_rot,
            nonrigid_max_trans=5
        )
        rbs.append(rb)

    return rbs

def run_86(
        pdb_file,
        xls
):
    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.domain.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    structures["MTOR"] = (pdb_file, "A", "orange", "J")
    structures["RICTOR"] = (pdb_file, "B", "blue", "K")
    structures["MLST8"] = (pdb_file, "C", "green", "L")
    structures["MSIN1"] = (pdb_file, "D", "red", "M")
    structures["CRIM"] = (pdb_file, "E", "purple", "N")
    structures["RBD"] = (pdb_file, "F", "brown", "O")
    structures["MSIN1PH"] = (pdb_file, "G", "tan", "P")
    structures["AKT1PH"] = (pdb_file, "H", "salmon", "Q")
    structures["KINASE"] = (pdb_file, "I", "black", "R")

    clones = dict()
    mols = dict()
    for component in structures.keys():
        print(component)
        pdb_file, chain, color, clone_chain = structures[component]

        mol = st.create_molecule(
            name=component,
            sequence=seqs[component],
            chain_id=chain
        )
        mols[component] = mol

        if component == "KINASE":
            res_ranges = [("N",(1,83)),("C",(90,272))]
        else:
            res_ranges = [("",(1,len(mol.get_residues())))]
        for termini, res_range in res_ranges:
            atomic_residues = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=chain,
                soft_check=True,
                res_range=res_range
            )

            mol.add_representation(
                residues=atomic_residues,
                resolutions=[1,10],
                color=color
            )

        # flex = list(mol.get_non_atomic_residues())
        # mol.add_representation(
        #     flex,
        #     resolutions=[10],
        #     color=color,
        #     setup_particles_as_densities=False
        # )

        clone = mol.create_clone(
            chain_id=clone_chain
        )
        clones[component] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    for component in mols.keys():
        mol_rb = create_rigid_body(
            mol=mols[component],
            component=component,
            dof=dof
        )

        clone_rb = create_rigid_body(
            mol=clones[component],
            component=component,
            dof=dof
        )

    output_objects = []

    # SYMMETRY
    center = IMP.algebra.Vector3D([174.619,175.188,177.687])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    for component in clones.keys():
        print(component)
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    m.update()

    # XLS
    for xl in xls:
        print(xl)
        if xl == "":
            continue
        # file, weight, type = xls[i]
        xl_r = get_xl_restraint(
            hier=root_hier,
            xl_file=xl+".csv",
            w=1,
            type="DSS"
        )
        xl_r.add_to_model()
        output_objects.append(xl_r)

    n_frames=1
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str("./output_{}".format(pdb_file.stem)),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=1,
        number_of_frames=n_frames,
        replica_exchange_maximum_temperature=5,
        score_moved=False
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    run_86(
        pdb_file=Path(Path.home(), "mtorc2/data/em/score_test/kinase_1.pdb"),
        xls=["all.domain"]
    )