from pathlib import Path
import time
import sys
import math
import pandas as pd
import argparse

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

sys.path.append(str(Path(Path.home(), "mtorc2/src")))
import params


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--param_file", type=str, required=True)
    parser.add_argument("--flex", action="store_true")
    parser.add_argument("--sym", action="store_true")
    parser.add_argument("--akt_tail", action="store_true")
    parser.add_argument("--phos_res", type=int)
    parser.add_argument("--xls", type=str)
    parser.add_argument("--map", type=str)
    parser.add_argument("--em_w", type=int)
    parser.add_argument("--ev", action="store_true")
    parser.add_argument("--conn", action="store_true")
    parser.add_argument("--shuffle", action="store_true")
    parser.add_argument("--n_frames", type=int, required=True)
    args = parser.parse_args()

    if args.phos_res and not args.akt_tail:
        raise RuntimeError("phos_res cannot be used without akt_tail")

    params.write_params(vars(args), Path(args.output_dir, "params.txt"))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    param_df = pd.read_csv(args.param_file)

    print(param_df.head())

    clones = dict()
    mols = dict()
    subunits = ["MTOR", "RICTOR", "MLST8", "MSIN1", "AKT1"]
    for subunit in subunits:
        print(subunit)
        subunit_row_ids = list(param_df[param_df["subunit"] == subunit].index)
        color = param_df.iloc[subunit_row_ids[0], param_df.columns.get_loc("color")]
        chain = param_df.iloc[subunit_row_ids[0], param_df.columns.get_loc("model_chain")]

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chain
        )
        mols[subunit] = mol

        print(subunit_row_ids)
        for row_id in subunit_row_ids:
            name = param_df.iloc[row_id]["name"]
            pdb_file = Path(param_df.iloc[row_id]["pdb_file"])
            pdb_chain = param_df.iloc[row_id]["pdb_chain"]
            start_pdb = param_df.iloc[row_id]["start_pdb"]
            end_pdb = param_df.iloc[row_id]["end_pdb"]
            offset = int(param_df.iloc[row_id]["offset"])
            em_prefix = param_df.iloc[row_id]["em_prefix"]
            res_per_comp = int(param_df.iloc[row_id]["res_per_comp"])
            print(name, pdb_file, pdb_chain, start_pdb, end_pdb, offset, em_prefix)

            atom_res = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=pdb_chain,
                soft_check=True,
                res_range=(start_pdb, end_pdb),
                offset=offset,
                model_num=1
            )

            if type(em_prefix) == str:
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=res_per_comp,
                    density_prefix=str(em_prefix),
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

        if args.flex:
            if subunit == "AKT1":
                print(mol.get_non_atomic_residues())
                flex_res = list(mol.get_non_atomic_residues())[:22]
            else:
                flex_res = list(mol.get_non_atomic_residues())

            mol.add_representation(
                flex_res,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

        if args.akt_tail and subunit == "AKT1":
            # Add any residues past >425 that aren't already represented.
            # Have to hard code these numbers in.
            akt_tail = mols["AKT1"].residue_range(
                a=str(426),
                b=str(480)
            )

            print(len(akt_tail))

            akt_tail_flex = [res for res in akt_tail if res not in mol.get_atomic_residues()]

            print(len(akt_tail_flex))

            mol.add_representation(
                akt_tail_flex,
                resolutions=[1],
                color=color
            )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        print(args.sym, type(args.sym))
        if args.sym:
            print("CREATING CLONE")
            clone = mol.create_clone(
                chain_id=chr(ord(chain) + 5)
            )
            clones[subunit] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    rbs = dict()
    for subunit in subunits:
        rbs[subunit] = list()
        copies = list()
        copies.append(mols[subunit])
        if args.sym:
            copies.append(clones[subunit])

        subunit_row_ids = list(param_df[param_df["subunit"] == subunit].index)
        for row_id in subunit_row_ids:
            rb_start = param_df.iloc[row_id]["rb_start"]
            rb_end = param_df.iloc[row_id]["rb_end"]
            rb_trans = param_df.iloc[row_id]["rb_trans"]
            rb_rot = param_df.iloc[row_id]["rb_rot"]
            non_rigid_trans = param_df.iloc[row_id]["flex_trans"]

            for copy in copies:
                rigid_parts = copy.residue_range(
                    a=str(rb_start),
                    b=str(rb_end)
                )

                non_rigid_parts = [res for res in copy.get_non_atomic_residues() if res in rigid_parts]

                rb_movers, rb = dof.create_rigid_body(
                    rigid_parts=rigid_parts,
                    nonrigid_parts=non_rigid_parts,
                    max_trans=rb_trans,
                    max_rot=rb_rot,
                    nonrigid_max_trans=non_rigid_trans
                )

                # if copy == mol:
                rbs[subunit].append(rb)

    # CONNECTIVITY
    output_objects = []
    if args.conn:
        for subunit in subunits:
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
                objects=mols[subunit],
                resolution=1
            )
            cr.add_to_model()
            output_objects.append(cr)

    print("args.akt_tail: {}".format(args.akt_tail))
    if args.akt_tail:
        print("AKT tail restraints")
        for i in range(425,480):
            conn_r = IMP.pmi.restraints.basic.DistanceRestraint(
                root_hier=root_hier,
                tuple_selection1=(i,i,"AKT1",0),
                tuple_selection2=(i+1,i+1,"AKT1",0),
                distancemin=3.5,
                distancemax=3.5,
                weight=100
            )
            conn_r.add_to_model()
            output_objects.append(conn_r)

    # EXCLUDED VOLUME
    ev_objects = list()
    for subunit in subunits:
        ev_objects.append(mols[subunit])

        if args.sym:
            ev_objects.append(clones[subunit])

    if args.ev:
        ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects=ev_objects,
            resolution=10
        )
        ev_r.add_to_model()
        output_objects.append(ev_r)

    # SYMMETRY
    center = IMP.algebra.Vector3D([166.9824, 166.9824, 166.9824])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    if args.sym:
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

    if args.xls:
        xls = args.xls.split(",")
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
    if args.map:
        sel = IMP.atom.Selection(
            hierarchy=root_hier,
            representation_type=IMP.atom.DENSITIES,
            copy_index=0
        )

        densities = sel.get_selected_particles()
        gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(
            densities,
            target_fn=str(args.map),
            scale_target_to_mass=True,
            # target_mass=25008.40380000003,
            slope=0.01,
            target_radii_scale=3.0,
            target_is_rigid_body=False,
            weight=args.em_w
        )
        gem.add_to_model()
        gem.set_label("BayesianEM")
        output_objects.append(gem)

    # ACTIVE SITE
    if args.phos_res:
        active_site_dr = IMP.pmi.restraints.basic.DistanceToPointRestraint(
            root_hier=s.get_hierarchy(),
            tuple_selection=(args.phos_res, args.phos_res, "AKT1", 0),
            anchor_point=IMP.algebra.Vector3D(160.938, 122.210, 187.189),
            radius=0,
            kappa=1
        )
        active_site_dr.add_to_model()
        output_objects.append(active_site_dr)

    # SHUFFLE
    if args.shuffle:
        excluded_rbs = list()
        excluded_rbs.extend(rbs["MTOR"])
        excluded_rbs.extend(rbs["RICTOR"])
        excluded_rbs.extend(rbs["MLST8"])
        excluded_rbs.extend(rbs["MSIN1"][:2])

        IMP.pmi.tools.shuffle_configuration(
            objects=root_hier,
            excluded_rigid_bodies=excluded_rbs,
            max_translation=100,
            max_rotation=2*math.pi,
            verbose=False,
            cutoff=10,
            niterations=25
        )

    print("Start optimization")
    dof.optimize_flexible_beads(100)

    print("Start sampling")
    rex = IMP.pmi.macros.ReplicaExchange(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(args.output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=10,
        number_of_frames=args.n_frames,
        replica_exchange_maximum_temperature=5
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/args.n_frames)