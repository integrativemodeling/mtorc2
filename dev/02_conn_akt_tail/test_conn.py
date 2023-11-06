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


def test_conn(
        output_dir,
        n_frames,
):
    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path("./3o96_align.pdb")
    structures["AKT1"] = list()
    # structures["AKT1"].append((pdb_file, "A", "purple", 1, 145, 1, 116, 0, None))
    structures["AKT1"].append((pdb_file, "A", "purple", 146, 480, 146, 429, 0, None))

    clones = dict()
    mols = dict()
    chain = "A"
    # for subunit in subunits.keys():
    for subunit in ["AKT1"]:
        print(subunit)
        # pdb_file, chain, color, clone_chain = structures[component]

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chain
        )
        mols[subunit] = mol

        for file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix in structures[subunit]:
            print(file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix)
            atom_res = mol.add_structure(
                pdb_fn=str(file),
                chain_id=struct_chain,
                soft_check=True,
                res_range=(start_struct, end_struct),
                offset=offset,
                model_num=1
            )
            print(atom_res)

            mol.add_representation(
                residues=atom_res,
                resolutions=[1,10],
                color=color
            )
        print(mol.get_non_atomic_residues())

        # mol.add_representation(
        #     list(mol.get_non_atomic_residues()),
        #     resolutions=[1],
        #     color=color,
        #     setup_particles_as_densities=False
        # )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        chain = chr(ord(chain) + 1)

    akt_tail = mols["AKT1"].residue_range(
        a=str(430),
        b=str(480)
    )
    mol.add_representation(
        akt_tail,
        resolutions=[1],
        color=color
    )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rb_ranges, rbs = dict(), dict()
    rb_ranges["AKT1"] = [(146, 480, 5)]

    for subunit in ["AKT1"]:
        rbs[subunit] = list()
        mol = mols[subunit]
        for rb_start, rb_end, rb_trans in rb_ranges[subunit]:
            print(subunit, rb_start, rb_end)
            rigid_parts = mol.residue_range(
                a=str(rb_start),
                b=str(rb_end)
            )

            rb_movers, rb = dof.create_rigid_body(
                rigid_parts=rigid_parts,
                nonrigid_parts=akt_tail,
                max_trans=0,
                max_rot=0,
                nonrigid_max_trans=5
            )
            rbs[subunit].append(rb)

    # dof.create_flexible_beads(
    #     flex_parts=akt_tail
    # )

    # CONNECTIVITY
    output_objects = []
    for i in range(429,480):
        conn_r = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=root_hier,
            tuple_selection1=(i,i,"AKT1",0),
            tuple_selection2=(i+1,i+1,"AKT1",0),
            distancemin=3.5,
            distancemax=3.5
        )
        conn_r.add_to_model()
        output_objects.append(conn_r)

    # conn_part = mol.residue_range(
    #     a=str(429),
    #     b=str(480)
    # )
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
        objects=mols["AKT1"][:428],
        disorderedlength=False,
        upperharmonic=False,
        resolution=1,
    )
    cr.add_to_model()
    output_objects.append(cr)

    # EXCLUDED VOLUME
    ev_objects = list()
    ev_objects.append(mols[subunit])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # ACTIVE SITE
    active_site_dr = IMP.pmi.restraints.basic.DistanceToPointRestraint(
        root_hier=s.get_hierarchy(),
        tuple_selection=(473, 473, "AKT1", 0),
        anchor_point=IMP.algebra.Vector3D(160.938, 122.210, 187.189),
        radius=0,
        kappa=1
    )
    active_site_dr.add_to_model()
    output_objects.append(active_site_dr)

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
        replica_exchange_maximum_temperature=5
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    test_conn(
        output_dir=Path("./output_0"),
        n_frames=1000
    )