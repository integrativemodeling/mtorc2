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


def sample_122(
        output_dir,
        em_comps,
        res_per_comp,
        flex,
        xls,
        map,
        em_w,
        shuffle,
        n_frames
):
    print("output_dir:      {}".format(output_dir))
    print("em_comps:        {}".format(em_comps))
    print("res_per_comp:    {}".format(res_per_comp))
    print("xls:             {}".format(xls))
    print("map:             {}".format(map))
    print("n_frames:        {}".format(n_frames))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    # subunits = {"MTOR": "A", "RICTOR": "B", "MLST8": "C", "MSIN1": "D", "AKT1": "E"}
    subunits = ["MTOR", "RICTOR", "MLST8", "MSIN1", "AKT1"]
    structures = dict()

    pdb_models_dir = Path(glob_data_dir, "em/models/J252")
    pdb_dir = Path(glob_data_dir, "pdb")

    structures["MTOR"] = [(Path(pdb_models_dir, "mtorc2_J252.pdb"), "A", "blue", 1, 2549, 1, 2549, 0, 0, "MTOR")]
    structures["RICTOR"] = [(Path(pdb_models_dir, "mtorc2_J252.pdb"), "E", "green", 1, 1708, 1, 1708, 0, 0, "RICTOR")]
    structures["MLST8"] = [(Path(pdb_models_dir, "mtorc2_J252.pdb"), "C", "white", 1, 326, 1, 326, 0, 0, "MLST8")]
    # [(1, 157), (158, 279), (280, 380), (381,522)]
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_models_dir, "mtorc2_J252.pdb"), "G", "red", 1, 155, 1, 155, 0, 0, "MSIN1"))
    structures["MSIN1"].append((Path(pdb_models_dir, "mSIN1_CRIM_AF-Q9BPZ7-F1_conf02.pdb"), "A", "red", 156, 279, 156, 267, 0, 5, "CRIM"))
    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 280, 380, 280, 353, 0, 5, "RBD"))
    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 381, 522, 381, 481, 0, 5, "MSIN1PH"))
    structures["AKT1"] = list()
    structures["AKT1"].append((Path(pdb_dir, "akt1_af.pdb"), "H", "purple", 1, 145, 1, 123, 0, 5, "AKT1PH"))
    structures["AKT1"].append((Path(pdb_models_dir, "AKT1_N-lobe_AF-P31749-F1_conf02.pdb"), "A", "purple", 146, 229, 146, 229, 0, 5, "KINASEN"))
    structures["AKT1"].append((Path(pdb_models_dir, "AKT1_C-lobe_AF-P31749-F1_conf02.pdb"), "A", "purple", 230, 480, 230, 410, 0, 5, "KINASEC"))

    clones = dict()
    mols = dict()
    chain = "A"
    # for subunit in subunits.keys():
    for subunit in subunits:
        print(subunit)
        # pdb_file, chain, color, clone_chain = structures[component]

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chain
        )
        mols[subunit] = mol

        for file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix in structures[subunit]:
            print(file, struct_chain, color)
            atom_res = mol.add_structure(
                pdb_fn=str(file),
                chain_id=struct_chain,
                soft_check=True,
                res_range=(start_struct, end_struct),
                offset=offset,
                model_num=1
            )

            if prefix in em_comps:
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=res_per_comp,
                    density_prefix=prefix,
                    density_force_compute=False,
                    resolutions=[1,10],
                    color=color
                )
            else:
                mol.add_representation(
                    residues=atom_res,
                    resolutions=[1,10],
                    color=color
                )

        # Don't add representation to the Akt tail starting at residue 411
        if flex and subunit == "AKT1":
            akt_main = list(mol.get_non_atomic_residues())[:22]
            akt_tail = list(mol.get_non_atomic_residues())[22:]
            akt_include = akt_main
            akt_include.extend(akt_tail[:16])

            print("AKT_MAIN: {}".format(akt_main))
            print("AKT_TAIL: {}".format(akt_tail))

            mol.add_representation(
                akt_include,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )
        elif flex:
            mol.add_representation(
                mol.get_non_atomic_residues(),
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        clone = mol.create_clone(
            chain_id=chr(ord(chain) + 5)
        )
        clones[subunit] = clone

        chain = chr(ord(chain) + 1)

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rb_ranges, rbs = dict(), dict()
    # rb_ranges["MTOR"] = [(1, len(mols["MTOR"].get_residues()), 0)]
    # rb_ranges["RICTOR"] = [(1, len(mols["RICTOR"].get_residues()), 0)]
    # rb_ranges["MLST8"] = [(1, len(mols["MLST8"].get_residues()), 0)]
    # rb_ranges["MSIN1"] = [(1, 155, 0), (156, 279, 5), (280, 380, 5), (381,522, 5)]
    # # rb_ranges["AKT1"] = [(1,144, 5), (145, 234, 5), (235, 480, 5)]
    # rb_ranges["AKT1"] = [(1, 145, 5), (146, 229, 5), (230, 480, 5)]

    for subunit in subunits:
        rbs[subunit] = list()
        mol, clone = mols[subunit], clones[subunit]

        for file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix in structures[subunit]:
        # for rb_start, rb_end, rb_trans in rb_ranges[subunit]:
            for copy in [mol, clone]:
                rigid_parts = copy.residue_range(
                    a=str(start),
                    b=str(end)
                )

                non_rigid_parts = [res for res in copy.get_non_atomic_residues() if res in rigid_parts]

                if copy == mol:
                    print(subunit, start, end)
                    print(rigid_parts)
                    print(non_rigid_parts)

                rb_movers, rb = dof.create_rigid_body(
                    rigid_parts=rigid_parts,
                    nonrigid_parts=non_rigid_parts,
                    max_trans=rb_trans,
                    max_rot=rb_trans,
                    nonrigid_max_trans=5
                )

                # if copy == mol:
                rbs[subunit].append(rb)

    # CONNECTIVITY
    output_objects = []
    for subunit in subunits:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[subunit],
            resolution=1
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for subunit in subunits:
        ev_objects.append(mols[subunit])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # SYMMETRY
    center = IMP.algebra.Vector3D([166.415,166.428,171.871])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    for subunit in subunits:
        dof.constrain_symmetry(
            references=mols[subunit],
            clones=clones[subunit],
            transform=transform
        )

    m.update()

    # XLS
    xl_dir = Path(Path.home(), "mtorc2/data/xlms/csvs")
    xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkwc.set_protein1_key("prot1")
    xldbkwc.set_protein2_key("prot2")
    xldbkwc.set_residue1_key("res1")
    xldbkwc.set_residue2_key("res2")

    xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)

    for xl in xls:
        xl, xl_w, xl_type = xl.split(":")
        xl_file = Path(xl_dir, xl+".csv")
        xldb.create_set_from_file(str(xl_file))

        if xl_type == "DSS":
            linker = ihm.cross_linkers.dss
            length = 35
        elif xl_type == "EDC":
            linker = ihm.cross_linkers.edc
            length = 16

        xl_r = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb,
            length=length,
            slope=0.01,
            resolution=1.,
            label=type,
            linker=linker,
            weight=float(xl_w)
        )

        xl_r.add_to_model()
        output_objects.append(xl_r)

    # EM
    sel = IMP.atom.Selection(
        hierarchy=root_hier,
        representation_type=IMP.atom.DENSITIES,
        copy_index=0
    )
    densities = sel.get_selected_particles()
    target_gmm_file = Path(glob_data_dir, "em/maps/{}".format(map))
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
        densities,
        target_fn=str(target_gmm_file),
        scale_target_to_mass=True,
        slope=0.01,
        target_radii_scale=3.0,
        target_is_rigid_body=False,
        weight=em_w
    )
    gem.add_to_model()
    gem.set_label("BayesianEM")
    output_objects.append(gem)

    # SHUFFLE
    excluded_rbs = list()
    excluded_rbs.extend(rbs["MTOR"])
    excluded_rbs.extend(rbs["RICTOR"])
    excluded_rbs.extend(rbs["MLST8"])
    excluded_rbs.extend(rbs["MSIN1"][:2])

    if shuffle:
        IMP.pmi.tools.shuffle_configuration(
            objects=root_hier,
            excluded_rigid_bodies=excluded_rbs,
            max_translation=100,
            max_rotation=2*math.pi,
            verbose=False,
            cutoff=10,
            niterations=25
        )
    dof.optimize_flexible_beads(100)

    rex = IMP.pmi.macros.ReplicaExchange(
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
    sample_122(
        output_dir=sys.argv[1],
        em_comps=sys.argv[2].split(","),
        res_per_comp=int(sys.argv[3]),
        flex=int(sys.argv[4]),
        xls=sys.argv[5].split(","),
        map=None if sys.argv[6] == "" else sys.argv[6],
        em_w=int(sys.argv[7]),
        shuffle=int(sys.argv[8]),
        n_frames=int(sys.argv[9])
    )
