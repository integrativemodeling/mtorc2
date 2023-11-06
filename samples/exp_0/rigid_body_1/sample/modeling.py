import IMP
import IMP.algebra
import IMP.pmi
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import ihm.cross_linkers
import math, time, os
import argparse
from pathlib import Path


if __name__ == "__main__":
    glob_data_dir = Path(Path.home(), "mtorc2/data")
    n_steps = 10000
    job_dir = Path(Path.home(), "mtorc2/rigid_body_1")

    fasta_file = Path(Path.home(), "mtorc2/data/fasta/MTOR2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    beadsize = 10

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)
    st = s.create_state()
    mols = dict()
    for component in ["RICTOR", "MTOR", "MLST8"]:
        mol = st.create_molecule(
            name=component,
            sequence=seqs[component]
        )
        mols[component] = mol
    print(mols)

    # Build Rictor
    # PDB files are 1-based indexed
    structures = dict()
    # Use RICTOR af structure because model and sequence are slightly off
    structures["RICTOR"] = (Path(glob_data_dir, "pdb/RICTOR.af.60.pdb"), "A", "blue")
    structures["MTOR"] = (Path(glob_data_dir, "pdb/MTOR.no_unk.pdb"), "A", "orange")
    structures["MLST8"] = (Path(glob_data_dir, "pdb/mTORC2.pdb"), "C", "green")
    clones = dict()
    for component in mols.keys():
        print(component)
        mol = mols[component]
        pdb_file, chain, color = structures[component]
        n_res = len(mol.get_residues())

        atom = mol.add_structure(
            pdb_fn=str(pdb_file),
            chain_id=chain
        )

        mol.add_representation(
             residues=atom,
             resolutions=[1,10],
             color=color
        )

        rictor_flex = mol.get_non_atomic_residues()
        mol.add_representation(
            rictor_flex,
            resolutions=[10],
            color=color
        )

        clone = mol.create_clone("B")
        clones[component] = clone

    root_hier = s.build()

    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    trans = IMP.algebra.Transformation3D([50,0,0])
    for component in mols.keys():
        mol = mols[component]
        clone = clones[component]
        # Small movements
        mov, rb = dof.create_rigid_body(
            mol,
            nonrigid_parts=mol.get_non_atomic_residues(),
            max_trans=1,
            max_rot=1,
            nonrigid_max_trans=1
        )

        mov_clone, rb_clone = dof.create_rigid_body(
            clone,
            nonrigid_parts=clone.get_non_atomic_residues(),
            max_trans=1,
            max_rot=1,
            nonrigid_max_trans=1
        )
        # IMP.atom.transform(rb_clone, trans)

        IMP.pmi.tools.display_bonds(mol)

    output_objects = []  # keep a list of functions that need to be reported
    crs = []

    for component in mols.keys():
        mol = mols[component]
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)

    evr_1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[mols[component] for component in mols.keys()],
        resolution=10
    )
    evr_1.add_to_model()
    output_objects.append(evr_1)

    evr_2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[mols["RICTOR"], clones["MTOR"]],
        resolution=10
    )
    evr_2.add_to_model()
    output_objects.append(evr_2)

    evr_3 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[mols["MTOR"], clones["MTOR"]],
        resolution=10
    )
    evr_3.add_to_model()
    output_objects.append(evr_3)

    center = IMP.algebra.Vector3D([196.362, 167.918, 128.446])
    for component in mols.keys():
        mol = mols[component]
        clone = clones[component]

        rot = IMP.algebra.get_rotation_about_axis([0, 0, 1], math.pi)
        transform = IMP.algebra.get_rotation_about_point(center, rot)
        dof.constrain_symmetry(mol, clone, transform)

    m.update()

    xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkwc.set_protein1_key("prot1")
    xldbkwc.set_protein2_key("prot2")
    xldbkwc.set_residue1_key("res1")
    xldbkwc.set_residue2_key("res2")

    xl_inter_file = Path(glob_data_dir, "xlms/xl_inter_formatted.csv")
    xldb_inter = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb_inter.create_set_from_file(str(xl_inter_file))

    xl_inter_r = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,
        database=xldb_inter,
        length=21.0,
        slope=0.01,
        resolution=1.0,
        label="INTER",
        linker=ihm.cross_linkers.dss,
        weight=1.
    )
    xl_inter_r.add_to_model()
    output_objects.append(xl_inter_r)

    xl_mlst8_mtor_file = Path(glob_data_dir, "xlms/xl_mlst8_mtor_formatted.csv")
    xl_mlst8_mtor_db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xl_mlst8_mtor_db.create_set_from_file(str(xl_mlst8_mtor_file))
    xl_mlst8_mtor_r = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,
        database=xl_mlst8_mtor_db,
        length=21.0,
        slope=0.01,
        resolution=1.0,
        label="INTER",
        linker=ihm.cross_linkers.dss,
        weight=2.
    )
    xl_mlst8_mtor_r.add_to_model()
    output_objects.append(xl_mlst8_mtor_r)

    # Don't shuffle components
    # IMP.pmi.tools.shuffle_configuration(root_hier, max_translation=30)
    # dof.optimize_flexible_beads(100)
    t0 = time.time()
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str("./output"),
        output_objects=output_objects,
        monte_carlo_steps=10,
        number_of_best_scoring_models=1,
        number_of_frames=n_steps)
    rex.execute_macro()
    print("TIME: {}".format(time.time() - t0))
