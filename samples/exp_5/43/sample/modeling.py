import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
from pathlib import Path
import sys, os
sys.path.append("/home/matthew/mtorc2/single_traj/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/single_traj/src")
from restraints import get_xl_restraint
import IMP.algebra
import math
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic


if __name__ == "__main__":
    n_frames = 100000

    if len(sys.argv) > 1:
        home_dir = Path(str(sys.argv[1]))
    else:
        home_dir = Path.home()

    exp_dir = Path(home_dir, "mtorc2/single_traj", Path.cwd().parents[1].name)
    output_dir = Path(exp_dir, Path.cwd().parents[0].name, "sample/output_0")
    print("OUT DIR:", output_dir)

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.rbd_ph.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(glob_data_dir, "pdb/em_dimer.kinase_crim.pdb")
    structures["MTOR"] = (pdb_file, "A", "orange", "J")
    structures["RICTOR"] = (pdb_file, "B", "blue", "K")
    structures["MLST8"] = (pdb_file, "C", "green", "L")
    structures["MSIN1"] = (pdb_file, "D", "red", "M")
    structures["CRIM"] = (pdb_file, "E", "purple", "N")
    structures["RBD"] = (Path(glob_data_dir, "pdb/msin1.rbd_ph_7lc1.pdb"), "F", "brown", "O")
    structures["AKT1"] = (Path(glob_data_dir, "pdb/AKT1.af.pdb"), "H", "salmon", "Q")

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

        # AKT1 has no flexible regions.
        if component not in ["AKT1"]:
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
    # RIGID BODIES
    for component in mols.keys():
        max_trans, max_rot = 1, 1
        if component not in ["AKT1"]:
            if component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
                max_trans, max_rot = 0, 0

            dof.create_rigid_body(
                rigid_parts=mols[component],
                nonrigid_parts=mols[component].get_non_atomic_residues(),
                max_trans=max_trans,
                max_rot=max_rot,
                nonrigid_max_trans=5
            )

            dof.create_rigid_body(
                rigid_parts=clones[component],
                nonrigid_parts=clones[component].get_non_atomic_residues(),
                max_trans=max_trans,
                max_rot=max_rot,
                nonrigid_max_trans=1
            )
        else:
            # AKT1 rigid body.
            dof.create_rigid_body(
                rigid_parts=mols["AKT1"],
                max_trans=5,
                max_rot=5
            )
            dof.create_rigid_body(
                rigid_parts=clones["AKT1"],
                max_trans=5,
                max_rot=5
            )

    output_objects = []
    # CONNECTIVITY
    connections = list()
    connections.append(("MSIN1", 157, "CRIM", 1))
    connections.append(("CRIM", 121, "RBD", 1))
    for connection in connections:
        prot1, res1, prot2, res2 = connection
        dr = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=s.get_hierarchy(),
            tuple_selection1=(res1, res1, prot1, 0),
            tuple_selection2=(res2, res2, prot2, 0),
            distancemin=0,
            distancemax=25,
            resolution=10,
            weight=1
        )
        dr.add_to_model()
        output_objects.append(dr)

    for component in mols.keys():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=10
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for component in mols.keys():
        ev_objects.append(mols[component])
        ev_objects.append(clones[component])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # SYMMETRY
    center = IMP.algebra.Vector3D([174.619,175.188,177.687])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    for component in ["CRIM", "RBD", "AKT1"]:
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    for component in ["MTOR", "RICTOR", "MSIN1", "MLST8"]:
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    m.update()

    # XLS
    xls = list()
    xls.append(("39.csv", 5, "DSS"))
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

    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=1,
        number_of_frames=n_frames
    )

    t0 = time.time()
    rex.execute_macro()
    print(time.time() - t0)

    if len(sys.argv) > 1:
        print("cp -r {} .".format(output_dir))
        os.system("cp -r {} .".format(output_dir))

