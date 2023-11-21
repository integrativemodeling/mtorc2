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
# sys.path.append("/wynton/home/sali/mhancock/mtorc2/src")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--param_file", type=str, required=True)
    parser.add_argument("--em_comps", type=str, required=True)
    parser.add_argument("--res_per_comp", type=int, required=True)
    parser.add_argument("--flex", type=int, required=True)
    parser.add_argument("--sym", type=str, required=True)
    parser.add_argument("--xls", type=str, required=True)
    parser.add_argument("--map", type=str, required=True)
    parser.add_argument("--em_w", type=int, required=True)
    parser.add_argument("--shuffle", type=int, required=True)
    parser.add_argument("--n_frames", type=int, required=True)
    args = parser.parse_args()

    params.write_params(vars(args), Path(args.output_dir, "params.txt"))

    em_comps = args.em_comps.split(",")

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
    subunits = list(set(param_df["subunit"]))
    # subunits = ["MTOR", "RICTOR", "MLST8", "MSIN1", "AKT1"]
    for subunit in subunits:
        print(subunit)
        # pdb_file, chain, color, clone_chain = structures[component]

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
            name = param_df.iloc[row_id, param_df.columns.get_loc("name")]
            pdb_file = param_df.iloc[row_id, param_df.columns.get_loc("pdb_file")]
            pdb_file = Path(glob_data_dir, pdb_file)
            pdb_chain = param_df.iloc[row_id, param_df.columns.get_loc("pdb_chain")]
            start_pdb = param_df.iloc[row_id, param_df.columns.get_loc("start_pdb")]
            end_pdb = param_df.iloc[row_id, param_df.columns.get_loc("end_pdb")]
            offset = param_df.iloc[row_id, param_df.columns.get_loc("offset")]

            print(name, pdb_file, pdb_chain, start_pdb, end_pdb, offset)
            atom_res = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=pdb_chain,
                soft_check=True,
                res_range=(start_pdb, end_pdb),
                offset=offset,
                model_num=1
            )

            if name in em_comps:
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=int(args.res_per_comp),
                    density_prefix=name,
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
                print(flex_res)
            else:
                flex_res = list(mol.get_non_atomic_residues())
            mol.add_representation(
                flex_res,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        if args.sym:
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
            rb_start = param_df.iloc[row_id, param_df.columns.get_loc("rb_start")]
            rb_end = param_df.iloc[row_id, param_df.columns.get_loc("rb_end")]
            rb_trans = param_df.iloc[row_id, param_df.columns.get_loc("rb_trans")]
            rb_rot = param_df.iloc[row_id, param_df.columns.get_loc("rb_rot")]
            non_rigid_trans = 5

            for copy in copies:
                rigid_parts = copy.residue_range(
                    a=str(rb_start),
                    b=str(rb_end)
                )

                non_rigid_parts = [res for res in copy.get_non_atomic_residues() if res in rigid_parts]

                if copy == mol:
                    print(subunit, rb_start, rb_end)
                    print(rigid_parts)
                    print(non_rigid_parts)

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

        if args.sym:
            ev_objects.append(clones[subunit])

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
    xls = args.xls.split(",")

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
            slope=0.01,
            target_radii_scale=3.0,
            target_is_rigid_body=False,
            weight=args.em_w
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

    if args.shuffle:
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