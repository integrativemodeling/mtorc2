import IMP
import IMP.algebra
# import IMP.pmi
# import IMP.pmi.io.crosslink
# import IMP.pmi.topology
# import IMP.pmi.dof
# import IMP.pmi.macros
# import IMP.pmi.restraints
# import IMP.pmi.restraints.stereochemistry
# import IMP.pmi.restraints.crosslinking
import math
from pathlib import Path


def get_representation(
        m,
        s,
        dimer=False,
        flex=False
):
    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/MTOR2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    mols = dict()
    for component in ["RICTOR", "MTOR", "MLST8"]:
        mol = st.create_molecule(
            name=component,
            sequence=seqs[component]
        )
        mols[component] = mol
    print(mols)

    # PDB files are 1-based indexed
    structures = dict()
    # Use RICTOR af structure because model and sequence are slightly off
    structures["MTOR"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "A", "orange", "D")
    structures["RICTOR"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "B", "blue", "E")
    structures["MLST8"] = (Path(glob_data_dir, "pdb/MTORC2.rigid_body.pdb"), "C", "green", "F")
    clones = dict()
    for component in mols.keys():
        print(component)
        mol = mols[component]
        pdb_file, chain, color, clone_chain = structures[component]

        atom = mol.add_structure(
            pdb_fn=str(pdb_file),
            chain_id=chain
        )

        mol.add_representation(
            residues=atom,
            resolutions=[1,10],
            color=color
        )

        if flex:
            flex = mol.get_non_atomic_residues()
            mol.add_representation(
                flex,
                resolutions=[10],
                color=color
            )

        if dimer:
            clone = mol.create_clone(
                chain_id=clone_chain
            )
            clones[component] = clone

        # IMP.pmi.tools.display_bonds(mol)

    if dimer:
        return mols, clones
    else:
        return mols


def create_rigid_bodies(
        dof,
        mols,
        clones=None,
):
    for component in mols.keys():
        mol = mols[component]
        dof.create_rigid_body(
            mol,
            nonrigid_parts=mol.get_non_atomic_residues(),
            max_trans=1,
            max_rot=1,
            nonrigid_max_trans=1
        )

        if clones:
            clone = clones[component]
            dof.create_rigid_body(
                clone,
                nonrigid_parts=clone.get_non_atomic_residues(),
                max_trans=1,
                max_rot=1,
                nonrigid_max_trans=1
            )


def setup_c2_symmetry(
        mols,
        clones,
        dof
):
    center = IMP.algebra.Vector3D([196.362, 167.918, 128.446])
    for component in mols.keys():
        mol = mols[component]
        clone = clones[component]

        rot = IMP.algebra.get_rotation_about_axis([0, 0, 1], math.pi)
        transform = IMP.algebra.get_rotation_about_point(center, rot)
        dof.constrain_symmetry(mol, clone, transform)