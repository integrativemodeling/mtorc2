import IMP.display
import IMP.pmi
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import ihm.cross_linkers
import sys, time
from pathlib import Path
import random
random.seed(1)


if __name__ == "__main__":
    job_id = 0
    n_steps = 10000

    home_dir = Path(Path.home())
    # home_dir = Path(Path.home(), "Documents")
    sys.path.append(str(Path(home_dir, "mtorc2/utility/misc")))
    from rb_functions import get_flex_regions

    fasta_file = Path(home_dir, "mtorc2/data/fasta/MTOR2.fasta")
    output_dir = Path(home_dir, "mtorc2/rigid_body.af/data/output")
    xl_file = Path(home_dir, "mtorc2/data/xlms/all_xls.csv")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    components = ["MTOR", "RICTOR", "MSIN1", "MLST8", "AKT1"]
    files = dict()
    pdb_path = Path(home_dir, "mtorc2/data/pdb")
    files["RICTOR"] = Path(pdb_path, "RICTOR.af.pdb")
    files["MSIN1"] = Path(pdb_path, "mSIN1.af.pdb")
    files["AKT1"] = Path(pdb_path, "AKT1.af.pdb")
    # files["MTOR"] = Path(pdb_path, "6ZW0.edit.pdb")
    files["MTOR"] = Path(pdb_path, "MTORC2.rigid_body.pdb")
    files["MLST8"] = Path(pdb_path, "MTORC2.rigid_body.pdb")
    chains = dict()
    chains["MTOR"] = "A"
    chains["RICTOR"] = "B"
    chains["MLST8"] = "C"
    chains["MSIN1"] = "G"
    chains["AKT1"] = "M"
    colors=dict()
    colors["MTOR"] = "orange"
    colors["RICTOR"] = "blue"
    colors["MLST8"] = "green"
    colors["MSIN1"] = "red"
    colors["AKT1"] = "white"

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)
    st = s.create_state()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    mols = dict()
    for component in components:
        mol = st.create_molecule(
            name=component,
            sequence=seqs[component],
            chain_id=chains[component]
        )
        mols[component] = mol


# Regions are left-inclusive/right-exclusive relative to the FASTA sequence (0-based indexing)

    rb_regions, flex_regions = dict(), dict()
    rb_regions["RICTOR"] = [[25, 502], [523, 637], [646, 857], [872, 1005], [1422, 1440], [1514, 1529], [1604, 1695]]
    rb_regions["MSIN1"] = [[0, 32], [84, 136], [136, 266], [278, 353]]
    rb_regions["AKT1"] =  [[1, 116], [145,445]]
    rb_regions["MTOR"] = [[0, len(mols["MTOR"].get_residues())]]
    rb_regions["MLST8"] = [[0, len(mols["MLST8"].get_residues())]]
    flex_regions["MLST8"] = [[0,7], [324,326]]
    flex_regions["MTOR"] = [[0,66], [74,81], [246,257], [289,385], [404,409], [466,477], [548,578], [633,643], [903,932], [1222,1260], [1814,1866], [2436,2491]]

    for component in components:
        mol = mols[component]
        # color = (random.random(), random.random(), random.random())
        rbs = rb_regions[component]
        for rb in rbs:
            start, end = rb[0]+1, rb[1]
            print(component)
            atom = mol.add_structure(
                pdb_fn=str(files[component]),
                res_range=[start,end],
                # chain_id=chains[component]
                chain_id=chains[component],
                soft_check=True
            )

            mol.add_representation(
                residues=atom,
                resolutions=[1,10],
                color=colors[component]
            )

        mol.add_representation(
            residues=mols[component].get_non_atomic_residues(),
            resolutions=[10],
            color=colors[component]
        )

    root_hier = s.build()

    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    for component in ["MTOR", "MLST8"]:
        mol = mols[component]
        dof.create_rigid_body(
            rigid_parts=mol,
            nonrigid_parts=mol.get_non_atomic_residues(),
            max_trans=3,
            max_rot=3,
            nonrigid_max_trans=5
        )

    for component in ["RICTOR", "MSIN1", "AKT1"]:
        mol = mols[component]
        n_res = len(mol.get_residues())
        flex_regions[component] = get_flex_regions(rb_regions[component], n_res)
        n_res = len(mol.get_residues())
        rbs = rb_regions[component]
        flexs = flex_regions[component]

        # Add rigid body region representation
        for rb in rbs:
            print(component, "RB", rb)
            start, end = rb[0], rb[1]
            # res_range is indexed by residue ID (1-based indexing)

            # print(atom)

            dof.create_rigid_body(
                # rigid_parts=mols["MTOR"],
                mol[start:end],
                # mol[rb[0]:rb[1]],
                max_trans=3,
                max_rot=3
                # nonrigid_max_trans=5
            )

        for flex in flexs:
            start, end = flex[0], flex[1]
            dof.create_flexible_beads(
                flex_parts=mol[start:end],
                max_trans=5,
                resolution=10
            )
            # del colors[0]

        # Add flexible region representation
        # for flex in flexs:
        #     print(component, "FLEX", flex)
        #     print(mol[flex[0]:flex[1]])
            # mol.add_representation(
            #     mol[flex[0]:flex[1]],
            #     resolutions=[10],
            #     color=colors[component]
            # )

    #
    # for component in ["RICTOR", "MSIN1", "AKT1"]:
    #     mol = mols[component]
    #     for rb in rb_regions[component]:
    #         dof.create_rigid_body(
    #             mol[rb[0]:rb[1]],
    #             max_trans=3,
    #             max_rot=3,
    #             nonrigid_max_trans=1
    #         )
    #     IMP.pmi.tools.display_bonds(mol)
    #
    # for component in ["MTOR", "MLST8"]:
    #     mol = mols[component]
    #     dof.create_rigid_body(
    #         mol,
    #         nonrigid_parts=mol.get_non_atomic_residues(),
    #         max_trans=3,
    #         max_rot=3,
    #         nonrigid_max_trans=1
    #     )
    # IMP.pmi.tools.display_bonds(mols["MTOR"])

    output_objects = []  # keep a list of functions that need to be reported
    crs = []
    mols_list = list()
    for component in components:
        mol = mols[component]
        mols_list.append(mol)

        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)

    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=mols_list,
        resolution=10
    )
    evr.add_to_model()
    output_objects.append(evr)

    xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkwc.set_protein1_key("prot1")
    xldbkwc.set_protein2_key("prot2")
    xldbkwc.set_residue1_key("res1")
    xldbkwc.set_residue2_key("res2")

    xl_db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xl_db.create_set_from_file(str(xl_file))

    intra_xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,
        database=xl_db,
        length=21.0,
        slope=0.01,
        resolution=1.0,
        label="INTRA",
        linker=ihm.cross_linkers.dss,
        weight=1.
    )

    intra_xl.add_to_model()
    output_objects.append(intra_xl)
    #
    # inter_xl_file = Path(xl_dir, "xl_inter_formatted.csv")
    # inter_xl_db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    # inter_xl_db.create_set_from_file(str(inter_xl_file))
    #
    # inter_xl = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    #     root_hier=root_hier,
    #     database=inter_xl_db,
    #     length=21.0,
    #     slope=0.01,
    #     resolution=1.0,
    #     label="INTER",
    #     linker=ihm.cross_linkers.dss,
    #     weight=1.
    # )
    #
    # inter_xl.add_to_model()
    # output_objects.append(inter_xl)

    t0 = time.time()
    IMP.pmi.tools.shuffle_configuration(root_hier, max_translation=30)
    # dof.optimize_flexible_beads(100)
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=25,
        number_of_best_scoring_models=1,
        number_of_frames=n_steps)
    rex.execute_macro()
    t1 = time.time()
    print(t1-t0)