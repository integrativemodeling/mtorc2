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


def model_akt_tail(
        output_dir,
        pdb_file,
        xls,
        phos_res,
        shuffle,
        n_frames
):
    print("output_dir:      {}".format(output_dir))
    print("pdb_file:        {}".format(pdb_file))
    print("xls:             {}".format(xls))
    print("phos_res:        {}".format(phos_res))
    print("shuffle:         {}".format(shuffle))
    print("n_frames:        {}".format(n_frames))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    param_file = Path(glob_data_dir, "params/122.csv")
    param_df = pd.read_csv(param_file)

    print(param_df.head())

    mols = dict()
    subunits = ["MTOR", "RICTOR", "MLST8", "MSIN1", "AKT1"]

    akt_tail = None

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
            name = param_df.iloc[row_id, param_df.columns.get_loc("name")]
            # pdb_file = param_df.iloc[row_id, param_df.columns.get_loc("pdb_file")]
            # pdb_file = Path(Path.home(), "mtorc2/manuscript/submission_2/models/126_0_1/cluster.0/cluster_center_model.pdb")
            pdb_file = Path(glob_data_dir, pdb_file)
            # pdb_chain = param_df.iloc[row_id, param_df.columns.get_loc("pdb_chain")]
            pdb_chain = chain
            start_pdb = param_df.iloc[row_id, param_df.columns.get_loc("start_pdb")]
            end_pdb = param_df.iloc[row_id, param_df.columns.get_loc("end_pdb")]
            offset = param_df.iloc[row_id, param_df.columns.get_loc("offset")]

            atom_res = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=pdb_chain,
                soft_check=True,
                res_range=(start_pdb, end_pdb),
                offset=offset,
                model_num=1
            )

            mol.add_representation(
                residues=atom_res,
                resolutions=[1,10],
                color=color
            )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        if subunit == "AKT1":
            akt_tail = mols["AKT1"].residue_range(
                a=str(411),
                b=str(480)
            )

            mol.add_representation(
                akt_tail,
                resolutions=[1],
                color=color
            )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    rbs = dict()
    for subunit in subunits:
        rbs[subunit] = list()
        mol = mols[subunit]
        subunit_row_ids = list(param_df[param_df["subunit"] == subunit].index)

        for row_id in subunit_row_ids:
            rb_start = param_df.iloc[row_id, param_df.columns.get_loc("rb_start")]
            rb_end = param_df.iloc[row_id, param_df.columns.get_loc("rb_end")]
            rb_trans = param_df.iloc[row_id, param_df.columns.get_loc("rb_trans")]
            rb_rot = param_df.iloc[row_id, param_df.columns.get_loc("rb_rot")]

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
                    max_trans=0,
                    max_rot=0,
                    nonrigid_max_trans=0
                )

                # if copy == mol:
                rbs[subunit].append(rb)

    dof.create_flexible_beads(
        flex_parts=akt_tail,
        max_trans=3,
        resolution="all"
    )

    # CONNECTIVITY
    output_objects = list()
    for i in range(410,480):
        conn_r = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=root_hier,
            tuple_selection1=(i,i,"AKT1",0),
            tuple_selection2=(i+1,i+1,"AKT1",0),
            distancemin=3.5,
            distancemax=3.5
        )
        conn_r.add_to_model()
        output_objects.append(conn_r)

    m.update()

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

    # SHUFFLE
    excluded_rbs = list()
    excluded_rbs.extend(rbs["MTOR"])
    excluded_rbs.extend(rbs["RICTOR"])
    excluded_rbs.extend(rbs["MLST8"])
    excluded_rbs.extend(rbs["MSIN1"])
    excluded_rbs.extend(rbs["AKT1"])

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

    # ACTIVE SITE
    active_site_dr = IMP.pmi.restraints.basic.DistanceToPointRestraint(
        root_hier=s.get_hierarchy(),
        tuple_selection=(phos_res, phos_res, "AKT1", 0),
        anchor_point=IMP.algebra.Vector3D(160.938, 122.210, 187.189),
        radius=0,
        kappa=1
    )
    active_site_dr.add_to_model()
    output_objects.append(active_site_dr)

    rex = IMP.pmi.macros.ReplicaExchange(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=100,
        number_of_frames=n_frames,
        replica_exchange_maximum_temperature=5
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    model_akt_tail(
        output_dir=sys.argv[1],
        pdb_file=sys.argv[2],
        xls=sys.argv[3].split(","),
        phos_res=int(sys.argv[4]),
        shuffle=int(sys.argv[5]),
        n_frames=int(sys.argv[6])
    )