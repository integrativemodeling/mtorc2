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
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import ihm.cross_linkers
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


def sample_114(
        output_dir
):
    print("output_dir:      {}".format(output_dir))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    structures["RICTOR"] = (Path("./model.pdb"), "B", "green", 1, 1708, 1, 1708, 0, None)
    structures["AKT1"] = (Path("./model.pdb"), "E", "purple", 1, 480, 1, 480, 0, None)

    mols = dict()
    for subunit in ["RICTOR", "AKT1"]:
        print(subunit)
        file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix = structures[subunit]
        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=struct_chain
        )
        mols[subunit] = mol

        atom_res = mol.add_structure(
            pdb_fn=str(file),
            chain_id=struct_chain,
            soft_check=True
        )
        print(len(atom_res))

        mol.add_representation(
            residues=atom_res,
            resolutions=[1,10],
            color=color
        )

        mol.add_representation(
            mol.get_non_atomic_residues(),
            resolutions=[10],
            color=color
        )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rb_ranges, rbs = dict(), dict()

    for subunit in ["RICTOR", "AKT1"]:
        rbs[subunit] = list()
        mol = mols[subunit]

        rb_movers, rb = dof.create_rigid_body(
            rigid_parts=mol.get_residues(),
            nonrigid_parts=mol.get_non_atomic_residues(),
            max_trans=0,
            max_rot=0,
            nonrigid_max_trans=5
        )
        rbs[subunit].append(rb)

    # CONNECTIVITY
    output_objects = []
    for component in ["RICTOR", "AKT1"]:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=1
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_objects = list()
    # for component in mols.keys():
    for subunit in ["RICTOR", "AKT1"]:
        ev_objects.append(mols[subunit])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    conn_dr = IMP.pmi.restraints.basic.DistanceRestraint(
        root_hier=root_hier,
        tuple_selection1=(11,11,"RICTOR",0),
        tuple_selection2=(322,322,"AKT1",0),
        distancemin=0,
        distancemax=16
    )
    conn_dr.add_to_model()
    output_objects.append(conn_dr)

    n_frames = 1000
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=10,
        number_of_frames=n_frames,
        replica_exchange_maximum_temperature=5
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    sample_114(
        output_dir=Path("./output_0")
    )