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


def sample_arch(
        output_dir,
        map,
        n_frames
):
    print("output_dir:      {}".format(output_dir))
    print("map:             {}".format(map))
    print("n_frames:        {}".format(n_frames))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    param_file = Path(glob_data_dir, "params/130.csv")
    param_df = pd.read_csv(param_file)

    print(param_df.head())

    clones = dict()
    mols = dict()
    # subunits = list(set(param_df["subunit"]))
    subunits = ["MSIN1", "AKT1"]
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

            if name in ["KINASEN", "KINASEC", "CRIM"]:
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=10,
                    density_prefix=name,
                    density_force_compute=False,
                    resolutions=[1,10],
                    color=color
                )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    rbs = dict()
    for subunit in subunits:
        rbs[subunit] = list()
        mol = mols[subunit]

        subunit_rows = param_df[param_df["subunit"] == subunit]
        subunit_rows = subunit_rows[(subunit_rows["name"] == "KINASEN") & (subunit_rows["name"] == "KINASEC") & (subunit_rows["name"] == "CRIM")]
        subunit_row_ids = list(subunit_rows.index)

        for row_id in subunit_row_ids:
            rb_start = param_df.iloc[row_id, param_df.columns.get_loc("rb_start")]
            rb_end = param_df.iloc[row_id, param_df.columns.get_loc("rb_end")]
            rb_trans = param_df.iloc[row_id, param_df.columns.get_loc("rb_trans")]
            rb_rot = param_df.iloc[row_id, param_df.columns.get_loc("rb_rot")]
            non_rigid_trans = 5

            for copy in [mol]:
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

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    m.update()

    # EM
    if map:
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
            weight=1
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

    print("Start sampling")
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
    sample_arch(
        output_dir=sys.argv[1],
        map=sys.argv[2],
        n_frames=int(sys.argv[3])
    )
