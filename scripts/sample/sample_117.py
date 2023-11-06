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


def sample_117(
        output_dir,
        n_frames,
        xls,
        map,
        active_site,
        em_comps
):
    print("output_dir:      {}".format(output_dir))
    print("n_frames:        {}".format(n_frames))
    print("xls:             {}".format(xls))
    print("map              {}".format(map))
    print("active_site      {}".format(active_site))
    print("em_comps         {}".format(active_site))

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
    structures["MTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "A", "blue", 1, 2549, 1, 2549, 0, None)]
    structures["RICTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "E", "green", 1, 1708, 1, 1708, 0, None)]
    structures["MLST8"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "C", "white", 1, 326, 1, 326, 0, None)]
    # [(1, 157), (158, 279), (280, 380), (381,522)]
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_dir, "mtorc2_238.pdb"), "G", "red", 1, 157, 1, 157, 0, None))
    structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, 381, -116, "CRIM"))
    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 280, 380, 280, 353, 0, "RBD"))
    structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 381, 522, 381, 481, 0, "MSIN1PH"))
    structures["AKT1"] = list()
    structures["AKT1"].append((Path(pdb_dir, "akt1_af.pdb"), "H", "purple", 1, 144, 1, 123, 0, "AKT1PH"))
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 145, 234, 145, 228, 0, "KINASEN"))
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 235, 480, 235, 429, 0, "KINASEC"))

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
                flex_density = False
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=5,
                    density_prefix=prefix,
                    density_force_compute=False,
                    resolutions=[1,10],
                    color=color
                )
            else:
                flex_density = False
                mol.add_representation(
                    residues=atom_res,
                    resolutions=[1,10],
                    color=color
                )

        if subunit == "AKT1":
            akt10 = list(mol.get_non_atomic_residues())[:50]
            print("AKT10: {}".format(akt10))
            akt1 = list(mol.get_non_atomic_residues())[51:]
            print("AKT1: {}".format(akt1))

            mol.add_representation(
                akt10,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

            mol.add_representation(
                akt1,
                resolutions=[1],
                color=color,
                setup_particles_as_densities=False
            )

            # flex_regions = list()
            # flex_regions.append((124,144,10))
            # flex_regions.append((229,234,10))
            # flex_regions.append((430,480,1))
            # for start, end, res in flex_regions:
            #     res_range = mol.residue_range(
            #         a=str(start),
            #         b=str(end)
            #     )
            #     mol.add_representation(
            #         res_range,
            #         resolutions=[res],
            #         color=color,
            #         setup_particles_as_densities=False
            #     )
        else:
            mol.add_representation(
                mol.get_non_atomic_residues(),
                resolutions=[10],
                color=color,
                setup_particles_as_densities=flex_density
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
    rb_ranges["AKT1"] = [(1,144, 5), (145, 234, 5), (235, 480, 5)]

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

                if subunit == "AKT1":
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
    for component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=10
        )
        cr.add_to_model()
        output_objects.append(cr)

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

    # EXCLUDED VOLUME
    ev_objects = list()
    # for component in mols.keys():
    for subunit in subunits:
        ev_objects.append(mols[subunit])
        ev_objects.append(clones[subunit])

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
            length = 30
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
    print("EM: {}".format(sel.get_selected_particles()))
    densities = sel.get_selected_particles()
    target_gmm_file = Path(glob_data_dir, "em/maps/{}.txt".format(map))
    gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
        densities,
        target_fn=str(target_gmm_file),
        scale_target_to_mass=True,
        slope=0.01,
        target_radii_scale=3.0,
        target_is_rigid_body=False
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
            radius=0,
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
    sample_117(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
        xls=sys.argv[3].split(","),
        map=None if sys.argv[4] == "" else sys.argv[4],
        active_site=int(sys.argv[5]),
        em_comps=sys.argv[6].split(",")
    )