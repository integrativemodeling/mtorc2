import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
from pathlib import Path
import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from restraints import get_ev_restraints
from restraints import get_xl_restraint
from restraints import get_crs
import IMP.algebra
import math


if __name__ == "__main__":
    output_dir = Path("./output_0")

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/MTOR2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    # PDB files are 1-based indexed
    structures = dict()
    # Use RICTOR af structure because model and sequence are slightly off
    structures["MTOR"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "A", "orange", "D")
    structures["RICTOR"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "B", "blue", "E")
    structures["MLST8"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "C", "green", "F")
    structures["MSIN1"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "G", "red", "H")
    clones = dict()
    mols = dict()
    for component in structures.keys():
        print(component)
        # mol = mols[component]
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

        clone = mol.create_clone(
            chain_id=clone_chain
        )
        clones[component] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    non_rigid_parts = list()
    rigid_parts = list()
    for component in mols.keys():
        rigid_parts.append(mols[component])
        non_rigid_parts.extend(mols[component].get_non_atomic_residues())
    dof.create_rigid_body(
        rigid_parts=rigid_parts,
        # nonrigid_parts=non_rigid_parts,
        max_trans=1,
        max_rot=1,
        nonrigid_max_trans=1
    )

    non_rigid_parts = list()
    rigid_parts = list()
    for component in clones.keys():
        rigid_parts.append(clones[component])
        non_rigid_parts.extend(clones[component].get_non_atomic_residues())
    dof.create_rigid_body(
        rigid_parts=rigid_parts,
        # nonrigid_parts=non_rigid_parts,
        max_trans=1,
        max_rot=1,
        nonrigid_max_trans=1
    )

    output_objects = []  # keep a list of functions that need to be reported
    crs = get_crs(
        mols=mols
    )
    for cr in crs:
        cr.add_to_model()
        output_objects.append(cr)

    ev_rs = get_ev_restraints(
        mols=mols,
        clones=clones
    )
    for ev_r in ev_rs:
        ev_r.add_to_model()
        output_objects.append(ev_r)

    center = IMP.algebra.Vector3D([196.362, 167.918, 128.446])
    for component in ["MTOR"]:
        mol = mols[component]
        clone = clones[component]

        rot = IMP.algebra.get_rotation_about_axis([0, 0, 1], math.pi)
        transform = IMP.algebra.get_rotation_about_point(center, rot)
        dof.constrain_symmetry(mol, clone, transform)

    m.update()

    xls = list()
    xls.append(("inter_dss.csv", 5, "DSS"))
    xls.append(("inter_edc.csv", 5, "EDC"))
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
        monte_carlo_steps=10,
        number_of_best_scoring_models=1,
        number_of_frames=25000,
        replica_exchange_swap=False
    )

    rex.execute_macro()