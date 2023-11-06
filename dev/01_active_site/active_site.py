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


def active_site(
        pdb_file,
        output_dir,
        n_frames
):
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
    structures["MTOR"] = [(pdb_file, "A", "blue", 1, 2549, 1, 2549, 0, None)]
    structures["RICTOR"] = [(pdb_file, "B", "green", 1, 1708, 1, 1708, 0, None)]
    structures["MLST8"] = [(pdb_file, "C", "white", 1, 326, 1, 326, 0, None)]
    # [(1, 157), (158, 279), (280, 380), (381,522)]
    structures["MSIN1"] = list()
    structures["MSIN1"].append((pdb_file, "D", "red", 1, 481, 1, 481, 0, None))
    structures["AKT1"] = list()
    structures["AKT1"].append((pdb_file, "E", "purple", 1, 480, 1, 429, 0, None))

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

            mol.add_representation(
                atom_res,
                resolutions=[1,10],
                color=color
            )

        if subunit == "AKT1":
            res_range = mols["AKT1"].residue_range(
                a=str(430),
                b=str(480)
            )

            mol.add_representation(
                res_range,
                resolutions=[1],
                color=color,
                setup_particles_as_densities=False
            )

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
    rb_ranges["MSIN1"] = [(1, len(mols["MSIN1"].get_residues()), 5)]
    rb_ranges["AKT1"] = [(1, len(mols["AKT1"].get_residues()), 5)]

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
                    max_trans=0,
                    max_rot=0,
                    nonrigid_max_trans=5
                )
                rbs[subunit].append(rb)

    # CONNECTIVITY
    output_objects = list()
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
    active_site(
        pdb_file=Path("./119_6_2_2_no_beads.pdb"),
        output_dir="./output",
        n_frames=1000
    )