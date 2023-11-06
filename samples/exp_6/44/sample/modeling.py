import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
from pathlib import Path
import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from restraints import get_xl_restraint
import IMP.algebra
import math
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic


if __name__ == "__main__":
    output_dir = Path("./output_0")

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.domain.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    # PDB files are 1-based indexed
    structures = dict()
    pdb_file = Path(glob_data_dir, "pdb/em_dimer.kinase_crim.pdb")
    structures["MTOR"] = (pdb_file, "A", "orange", "J")
    structures["RICTOR"] = (pdb_file, "B", "blue", "K")
    structures["MLST8"] = (pdb_file, "C", "green", "L")
    structures["MSIN1"] = (pdb_file, "D", "red", "M")
    structures["CRIM"] = (pdb_file, "E", "purple", "N")
    # structures["RBD"] = (Path(glob_data_dir, "pdb/rbd.pdb"), "F", "brown", "O")
    # structures["MSIN1PH"] = (Path(glob_data_dir, "pdb/msin1.ph.pdb"), "G", "tan", "P")
    # structures["AKT1PH"] = (Path(glob_data_dir, "pdb/akt1.ph.pdb"), "H", "salmon", "Q")
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

        atom = mol.add_structure(
            pdb_fn=str(pdb_file),
            chain_id=chain,
            soft_check=True
        )

        mol.add_representation(
            residues=atom,
            resolutions=[1,10],
            color=color
        )
        flex = mol.get_non_atomic_residues()
        mol.add_representation(
            flex,
            resolutions=[10],
            color=color
        )

        clone = mol.create_clone(
            chain_id=clone_chain
        )
        clones[component] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # Main rigid body.
    for component in ["MTOR", "RICTOR", "MLST8", "MSIN1", "CRIM", "KINASE"]:
        dof.create_rigid_body(
            rigid_parts=mols[component],
            nonrigid_parts=mols[component].get_non_atomic_residues(),
            max_trans=0,
            max_rot=0,
            nonrigid_max_trans=1
        )

        dof.create_rigid_body(
            rigid_parts=clones[component],
            nonrigid_parts=clones[component].get_non_atomic_residues(),
            max_trans=0,
            max_rot=0,
            nonrigid_max_trans=1
        )

    # # CRIM, RBD rigid bodies.
    # for component in ["CRIM", "RBD"]:
    #     dof.create_rigid_body(
    #         rigid_parts=[mols[component]],
    #         nonrigid_parts=[mols[component].get_non_atomic_residues()],
    #         max_trans=5,
    #         max_rot=5,
    #         nonrigid_max_trans=1
    #     )
    #
    #     dof.create_rigid_body(
    #         rigid_parts=[clones[component]],
    #         nonrigid_parts=[clones[component].get_non_atomic_residues()],
    #         max_trans=5,
    #         max_rot=5,
    #         nonrigid_max_trans=1
    #     )
    #
    # # AKT1 rigid body.
    # dof.create_rigid_body(
    #     rigid_parts=mols["AKT1"],
    #     max_trans=5,
    #     max_rot=5
    # )
    # dof.create_rigid_body(
    #     rigid_parts=clones["AKT1"],
    #     max_trans=5,
    #     max_rot=5
    # )
    #
    output_objects = []  # keep a list of functions that need to be reported
    # # CONNECTIVITY
    # d_0 = IMP.pmi.restraints.basic.DistanceRestraint(
    #     root_hier=s.get_hierarchy(),
    #     tuple_selection1=(138, 138, "MSIN1"),
    #     tuple_selection2=(1, 1, "CRIM"),
    #     distancemax=10,
    #     resolution=10,
    #     weight=100
    # )
    # d_0.add_to_model()
    # output_objects.append(d_0)
    #
    # d_1 = IMP.pmi.restraints.basic.DistanceRestraint(
    #     root_hier=s.get_hierarchy(),
    #     tuple_selection1=(128, 128, "CRIM"),
    #     tuple_selection2=(1, 1, "RBD"),
    #     distancemax=10,
    #     resolution=10,
    #     weight=100
    # )
    # d_1.add_to_model()
    # output_objects.append(d_1)
    #
    for component in mols.keys():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=10
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_rs = list()

    r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[mols["MTOR"], mols["RICTOR"], mols["MLST8"], mols["MSIN1"], mols["CRIM"], mols["KINASE"], clones["MTOR"], clones["RICTOR"], clones["MLST8"], clones["MSIN1"], clones["CRIM"], clones["KINASE"]],
        resolution=10
    )

    ev_rs.append(r)
    #
    # for component in ["MSIN1", "CRIM", "RBD", "MTOR", "RICTOR", "MLST8"]:
    #     r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    #         included_objects=[mols["AKT1"], mols[component]],
    #         resolution=10
    #     )
    #     ev_rs.append(r)
    #
    #     r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    #         included_objects=[mols["AKT1"], clones[component]],
    #         resolution=10
    #     )
    #     ev_rs.append(r)
    #
    # r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    #     included_objects=[mols["AKT1"], clones["AKT1"]],
    #     resolution=10
    # )
    # ev_rs.append(r)

    for ev_r in ev_rs:
        ev_r.add_to_model()
        output_objects.append(ev_r)

    # SYMMETRY
    center = IMP.algebra.Vector3D([174.619,175.188,177.687])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    for component in ["MTOR", "RICTOR", "MSIN1", "MLST8", "CRIM", "KINASE"]:
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    m.update()

    # XLS
    xls = list()
    xls.append(("44.csv", 5, "DSS"))

    for i in range(len(xls)):
        file, weight, type = xls[i]
        xl_r = get_xl_restraint(
            hier=root_hier,
            xl_file=file,
            w=weight,
            type=type
        )
        xl_r.add_to_model()
        output_objects.append(xl_r)

    # Don't shuffle components
    t0 = time.time()
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=1,
        number_of_frames=10000
    )

    rex.execute_macro()
    print(time.time() - t0)
