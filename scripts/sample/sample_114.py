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
        output_dir,
        res_per_comp,
        akt_tail_res,
        em_comps,
        full_arch_ev,
        xls,
        map,
        em_w,
        active_site,
        n_frames
):
    print("output_dir:      {}".format(output_dir))
    print("res_per_comp:    {}".format(res_per_comp))
    print("akt_tail_res:    {}".format(akt_tail_res))
    print("em_comps:        {}".format(em_comps))
    print("full_arch_ev:    {}".format(full_arch_ev))
    print("xls:             {}".format(xls))
    print("map:             {}".format(map))
    print("em_w:            {}".format(em_w))
    print("active_site:     {}".format(active_site))
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
    pdb_dir = Path(Path.home(), "mtorc2/data/pdb")
    structures["MTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "A", "blue", 1, 2549, 1, 2549, 0, "MTOR")]
    structures["RICTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "E", "green", 1, 1708, 1, 1708, 0, "RICTOR")]
    structures["MLST8"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "C", "white", 1, 326, 1, 326, 0, "MLST8")]
    # [(1, 157), (158, 279), (280, 380), (381,522)]
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_dir, "mtorc2_238.pdb"), "G", "red", 1, 157, 1, 157, 0, "MSIN1"))
    # structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, 381, -116, "CRIM"))
    structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, 277+116, -116, "CRIM"))

    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 280, 380, 280, 353, 0, "RBD"))
    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 381, 522, 381, 481, 0, "MSIN1PH"))
    structures["AKT1"] = list()
    structures["AKT1"].append((Path(pdb_dir, "akt1_af.pdb"), "H", "purple", 1, 144, 1, 123, 0, "AKT1PH"))
    # structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 145, 234, 145, 228, 0, "KINASEN"))
    # structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 235, 480, 235, 429, 0, "KINASEC"))
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 145, 231, 145, 231, 0, "KINASEN"))
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 232, 480, 232, 429, 0, "KINASEC"))

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

        if subunit == "AKT1":
            akt_main = list(mol.get_non_atomic_residues())[:50]
            print("AKT_MAIN: {}".format(akt_main))
            akt_tail = list(mol.get_non_atomic_residues())[51:]
            print("AKT_TAIL: {}".format(akt_tail))

            mol.add_representation(
                akt_main,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

            mol.add_representation(
                akt_tail,
                resolutions=[akt_tail_res],
                color=color,
                setup_particles_as_densities=False
            )
        else:
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
    rb_ranges["MTOR"] = [(1, len(mols["MTOR"].get_residues()), 0)]
    rb_ranges["RICTOR"] = [(1, len(mols["RICTOR"].get_residues()), 0)]
    rb_ranges["MLST8"] = [(1, len(mols["MLST8"].get_residues()), 0)]
    rb_ranges["MSIN1"] = [(1, 157, 0), (158, 279, 5), (280, 380, 5), (381,522, 5)]
    # rb_ranges["AKT1"] = [(1,144, 5), (145, 234, 5), (235, 480, 5)]
    rb_ranges["AKT1"] = [(1,144, 5), (145, 231, 5), (232, 480, 5)]

    for subunit in subunits:
        rbs[subunit] = list()
        mol, clone = mols[subunit], clones[subunit]
        for rb_start, rb_end, rb_trans in rb_ranges[subunit]:
            for copy in [mol, clone]:
                print(subunit, rb_start, rb_end)
                rigid_parts = copy.residue_range(
                    a=str(rb_start),
                    b=str(rb_end)
                )

                non_rigid_parts = [res for res in copy.get_non_atomic_residues() if res in rigid_parts]

                print(rigid_parts)
                print(non_rigid_parts)

                rb_movers, rb = dof.create_rigid_body(
                    rigid_parts=rigid_parts,
                    nonrigid_parts=non_rigid_parts,
                    max_trans=rb_trans,
                    max_rot=rb_trans,
                    nonrigid_max_trans=5
                )
                rbs[subunit].append(rb)

    # CONNECTIVITY
    output_objects = []
    cr_main_comps = ["MTOR", "RICTOR", "MLST8", "MSIN1"]

    if akt_tail_res == 1:
        cr_akt = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols["AKT1"][:429],
            resolution=10
        )
        cr_akt.add_to_model()
        output_objects.append(cr_akt)

        for i in range(429,480):
            conn_dr = IMP.pmi.restraints.basic.DistanceRestraint(
                root_hier=root_hier,
                tuple_selection1=(i,i,"AKT1",0),
                tuple_selection2=(i+1,i+1,"AKT1",0),
                distancemin=3.5,
                distancemax=3.5
            )
            conn_dr.add_to_model()
            output_objects.append(conn_dr)
    else:
        cr_main_comps.append("AKT1")

    for component in cr_main_comps:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=1
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for subunit in subunits:
        for copy in mols, clones:
            if not full_arch_ev and subunit == "Akt1":
                ev_objects.append(
                    copy[subunit].residue_range(
                        a=str(1),
                        b=str(144)
                    )
                )

                ev_objects.append(
                    copy[subunit].residue_range(
                        a=str(232),
                        b=str(480)
                    )
                )
            else:
                ev_objects.append(
                    copy[subunit].get_residues()
                )

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects=ev_objects,
            resolution=10
        )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    if not full_arch_ev:
        arch_ev_partners = [("MTOR", 0), ("RICTOR",0), ("MLST8", 0), ("MSIN1", 0), ("MSIN1", 2), ("MSIN1", 3), ("AKT1", 0)]
        arch_ev_objects = list()
        for subunit, rb_num in arch_ev_partners:
            start, end, tmp = rb_ranges[subunit][rb_num]

            arch_ev_objects.append(
                mols["AKT1"].residue_range(
                    a=str(145),
                    b=str(231)
                )
            )

            arch_ev_objects.append(
                mols[subunit].residue_range(
                    a=str(start),
                    b=str(end)
                )
            )

            ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                included_objects=arch_ev_objects,
                resolution=10
            )
            ev_r.add_to_model()
            output_objects.append(ev_r)


# for subunit in subunits:
    #     for i in range(len(rbs[subunit])):
    #         if not full_arch_ev and subunit == "AKT1" and i == 1:
    #             continue
    #         else:
    #             print(rbs[subunit][i])
    #             ev_objects.append(rbs[subunit][i])
    #

    #
    # # Don't add Crim
    # if not full_arch_ev:
    #     ev_arch_objects = list()
    #     ev_arch_objects.append(rbs["AKT1"][1])
    #
    #     ev_arch_objects.append(rbs["MTOR"][0])
    #     ev_arch_objects.append(rbs["RICTOR"][0])
    #     ev_arch_objects.append(rbs["MLST8"][0])
    #
    #     # No Crim or Kinase C lobe
    #     ev_arch_objects.append(rbs["MSIN1"][0])
    #     ev_arch_objects.append(rbs["MSIN1"][2])
    #     ev_arch_objects.append(rbs["MSIN1"][3])
    #     ev_arch_objects.append(rbs["AKT1"][0])
    #
    #     arch_ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    #         included_objects=ev_arch_objects,
    #         resolution=10
    #     )
    #     arch_ev_r.add_to_model()
    #     output_objects.append(arch_ev_r)

    # SYMMETRY
    # center = IMP.algebra.Vector3D([174.619,175.188,177.687])
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
    xl_dir = Path(Path.home(), "mtorc2/data/xlms")
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
    target_gmm_file = Path(glob_data_dir, "em/maps/{}.txt".format(map))
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

    # ACTIVE SITE
    if active_site:
        active_site_dr = IMP.pmi.restraints.basic.DistanceToPointRestraint(
            root_hier=s.get_hierarchy(),
            tuple_selection=(473, 473, "AKT1", 0),
            anchor_point=IMP.algebra.Vector3D(160.938, 122.210, 187.189),
            radius=15,
            kappa=1
        )
        active_site_dr.add_to_model()
        output_objects.append(active_site_dr)

    # SHUFFLE
    shuffle_rbs = list()
    shuffle_rbs.extend(rbs["MSIN1"][1:])
    shuffle_rbs.extend(rbs["AKT1"])
    IMP.pmi.tools.shuffle_configuration(
        objects=shuffle_rbs,
        max_translation=100,
        max_rotation=2*math.pi,
        verbose=False,
        cutoff=10,
        niterations=25
    )
    dof.optimize_flexible_beads(100)

    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=10,
        number_of_frames=n_frames,
        replica_exchange_maximum_temperature=5,
        score_moved=False
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    sample_114(
        output_dir=sys.argv[1],
        res_per_comp=int(sys.argv[2]),
        akt_tail_res=int(sys.argv[3]),
        em_comps=sys.argv[4].split(","),
        full_arch_ev=int(sys.argv[5]),
        xls=sys.argv[6].split(","),
        map=None if sys.argv[7] == "" else sys.argv[7],
        em_w=int(sys.argv[8]),
        active_site=int(sys.argv[9]),
        n_frames=int(sys.argv[10])
    )