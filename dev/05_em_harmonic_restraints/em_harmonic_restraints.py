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


def em_harmonic_restraints(
        n_frames
):
    print("n_frames:        {}".format(n_frames))
    print("map              {}".format(map))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    # subunits = {"MTOR": "A", "RICTOR": "B", "MLST8": "C", "MSIN1": "D", "AKT1": "E"}
    subunits = ["MTOR", "RICTOR", "MLST8", "MSIN1"]
    structures = dict()
    pdb_dir = Path(Path.home(), "mtorc2/data/pdb")
    structures["MTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "A", "blue", 1, 2549, 1, 2549, 0, "MTOR")]
    structures["RICTOR"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "E", "green", 1, 1708, 1, 1708, 0, "RICTOR")]
    structures["MLST8"] = [(Path(pdb_dir, "mtorc2_238.pdb"), "C", "white", 1, 326, 1, 326, 0, "MLST8")]
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_dir, "mtorc2_238.pdb"), "G", "red", 1, 95, 1, 95, 0, "MSIN1"))
    structures["MSIN1"].append((Path(pdb_dir, "mtorc2_238.pdb"), "G", "red", 96, 157, 96, 157, 0, "MSIN1"))

    # structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, 381, -116, "CRIM"))
    # structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 280, 380, 280, 353, 0, "RBD"))
    # structures["MSIN1"].append((Path(pdb_dir, "7lc1.pdb"), "B", "red", 381, 522, 381, 481, 0, "MSIN1PH"))
    # structures["AKT1"] = list()
    # structures["AKT1"].append((Path(pdb_dir, "akt1_af.pdb"), "H", "purple", 1, 144, 1, 123, 0, "AKT1PH"))
    # structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 145, 234, 145, 228, 0, "KINASEN"))
    # structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 235, 480, 235, 429, 0, "KINASEC"))

    # clones = dict()
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
            print(len(atom_res))

            mol.add_representation(
                residues=atom_res,
                resolutions=[1,10],
                color=color
            )

        mol.add_representation(
            mol.get_non_atomic_residues(),
            resolutions=[10],
            color=color
        )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

        # clone = mol.create_clone(
        #     chain_id=chr(ord(chain) + 5)
        # )
        # clones[subunit] = clone

        chain = chr(ord(chain) + 1)

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rb_ranges, rbs = dict(), dict()
    rb_ranges["MTOR"] = [(1, len(mols["MTOR"].get_residues()), 0)]
    rb_ranges["RICTOR"] = [(1, len(mols["RICTOR"].get_residues()), 0)]
    rb_ranges["MLST8"] = [(1, len(mols["MLST8"].get_residues()), 0)]
    # rb_ranges["MSIN1"] = [(1, 157, 0)]

    # rb_ranges["MSIN1"] = [(1, 157, 0), (158, 279, 5), (280, 380, 5), (381,522, 5)]
    # rb_ranges["AKT1"] = [(1,144, 5), (145, 234, 5), (235, 480, 5)]

    msin1_rictor = list(mols["MSIN1"].residue_range(
        a=str(1),
        b=str(95)
    ))
    rictor_rb = list(mols["RICTOR"].get_residues())
    rictor_rb.extend(msin1_rictor)
    rb_movers, rb = dof.create_rigid_body(
        rigid_parts=rictor_rb,
        nonrigid_parts=mols["RICTOR"].get_non_atomic_residues(),
        max_trans=1,
        max_rot=1,
        nonrigid_max_trans=5
    )

    msin1_mlst8 = list(mols["MSIN1"].residue_range(
        a=str(96),
        b=str(157)
    ))
    mlst8_rb = list(mols["MLST8"].get_residues())
    mlst8_rb.extend(msin1_mlst8)

    print(mlst8_rb)

    rb_movers, rb = dof.create_rigid_body(
        rigid_parts=mlst8_rb,
        nonrigid_parts=mols["MLST8"].get_non_atomic_residues(),
        max_trans=1,
        max_rot=1,
        nonrigid_max_trans=5
    )

    mtor_rb = list(mols["MTOR"].get_residues())
    rb_movers, rb = dof.create_rigid_body(
        rigid_parts=mtor_rb,
        nonrigid_parts=mols["MTOR"].get_non_atomic_residues(),
        max_trans=1,
        max_rot=1,
        nonrigid_max_trans=5
    )

    # CONNECTIVITY
    output_objects = []
    for component in mols.keys():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=1
        )
        cr.add_to_model()
        output_objects.append(cr)

    # EXCLUDED VOLUME
    # ev_objects = list()
    # # for component in mols.keys():
    # for subunit in subunits:
    #     ev_objects.append(mols[subunit])
    #     # ev_objects.append(clones[subunit])

    ev_objects = list()
    # ev_objects.append(mols["MTOR"])
    # ev_objects.append(mols["RICTOR"])
    # ev_objects.append(mols["MLST8"])
    ev_objects.append(rictor_rb)
    ev_objects.append(mtor_rb)
    ev_objects.append(mlst8_rb)

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    m.update()

    # EM anchors
    em_anchor_1 = IMP.pmi.restraints.basic.DistanceToPointRestraint(
        root_hier=root_hier,
        tuple_selection=(436,436,"RICTOR",0),
        anchor_point=IMP.algebra.Vector3D(106.466,137.072,162.367),
        radius=0,
        kappa=1,
        resolution=1
    )
    em_anchor_1.add_to_model()
    output_objects.append(em_anchor_1)

    em_anchor_2 = IMP.pmi.restraints.basic.DistanceToPointRestraint(
        root_hier=root_hier,
        tuple_selection=(874,874,"MTOR",0),
        anchor_point=IMP.algebra.Vector3D(196.785,140.251,149.861),
        radius=0,
        kappa=1,
        resolution=1
    )
    em_anchor_2.add_to_model()
    output_objects.append(em_anchor_2)

    em_anchor_3 = IMP.pmi.restraints.basic.DistanceToPointRestraint(
        root_hier=root_hier,
        tuple_selection=(177,177,"MLST8",0),
        anchor_point=IMP.algebra.Vector3D(187.282,87.795,213.710),
        radius=0,
        kappa=1,
        resolution=1
    )
    em_anchor_3.add_to_model()
    output_objects.append(em_anchor_3)

    # EM restraints
    # em_dr_1_1 = IMP.pmi.restraints.basic.DistanceRestraint(
    #     root_hier=root_hier,
    #     tuple_selection1=(314,314,"MLST8",0),
    #     tuple_selection2=(2278,2278,"MTOR",0),
    #     distancemin=12.734819983,
    #     distancemax=12.734819983
    # )
    # em_dr_1_1.add_to_model()
    # output_objects.append(em_dr_1_1)

    def get_em_dr(
            chain_1,
            res_id_1,
            chain_2,
            res_id_2
    ):
        chain_names = {"A": "MTOR", "B": "RICTOR", "C": "MLST8"}
        sel_1 = IMP.atom.Selection(hierarchy=root_hier,chain_id=chain_1,residue_index=res_id_1)
        p_1 = sel_1.get_selected_particles()[0]
        xyz_1 = IMP.core.XYZ(p_1)

        sel_2 = IMP.atom.Selection(hierarchy=root_hier,chain_id=chain_2,residue_index=res_id_2)
        p_2 = sel_2.get_selected_particles()[0]
        xyz_2 = IMP.core.XYZ(m, p_2)
        d = IMP.algebra.get_distance(xyz_1.get_coordinates(), xyz_2.get_coordinates())
        dr = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=root_hier,
            tuple_selection1=(res_id_1,res_id_1,chain_names[chain_1],0),
            tuple_selection2=(res_id_2,res_id_2,chain_names[chain_2],0),
            distancemin=d,
            distancemax=d
        )

        return dr

    dr_1 = get_em_dr(chain_1="A", res_id_1=2515, chain_2="C", res_id_2=169)
    dr_1.add_to_model()
    output_objects.append(dr_1)

    dr_2 = get_em_dr(chain_1="A", res_id_1=2278, chain_2="C", res_id_2=314)
    dr_2.add_to_model()
    output_objects.append(dr_2)

    dr_3 = get_em_dr(chain_1="B", res_id_1=1569, chain_2="C", res_id_2=118)
    dr_3.add_to_model()
    output_objects.append(dr_3)

    dr_4 = get_em_dr(chain_1="B", res_id_1=1560, chain_2="C", res_id_2=84)
    dr_4.add_to_model()
    output_objects.append(dr_4)

    dr_5 = get_em_dr(chain_1="A", res_id_1=1217, chain_2="B", res_id_2=558)
    dr_5.add_to_model()
    output_objects.append(dr_3)

    dr_6 = get_em_dr(chain_1="A", res_id_1=2076, chain_2="B", res_id_2=248)
    dr_6.add_to_model()
    output_objects.append(dr_4)


    # em_dr_2 = IMP.pmi.restraints.basic.DistanceRestraint(
    #     root_hier=root_hier,
    #     tuple_selection1=(118,118,"MLST8",0),
    #     tuple_selection2=(1569,1569,"RICTOR",0),
    #     distancemin=13.2460217424,
    #     distancemax=13.2460217424
    # )
    # em_dr_2.add_to_model()
    # output_objects.append(em_dr_2)
    #
    # em_dr_3 = IMP.pmi.restraints.basic.DistanceRestraint(
    #     root_hier=root_hier,
    #     tuple_selection1=(1217,1217,"MTOR",0),
    #     tuple_selection2=(558,558,"RICTOR",0),
    #     distancemin=9.94589608834,
    #     distancemax=9.94589608834
    # )
    # em_dr_3.add_to_model()
    # output_objects.append(em_dr_3)

    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(Path("./output_0")),
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
    em_harmonic_restraints(
        n_frames=100
    )