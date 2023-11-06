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


def test_em_mtorc2(
        n_frames,
        map
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
                density_residues_per_component=10,
                density_prefix=prefix,
                density_force_compute=False,
                resolutions=[1,10],
                color=color
            )

        mol.add_representation(
            mol.get_non_atomic_residues(),
            resolutions=[10],
            color=color,
            setup_particles_as_densities=False
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
    ev_objects = list()
    # for component in mols.keys():
    for subunit in subunits:
        ev_objects.append(mols[subunit])
        # ev_objects.append(clones[subunit])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # # SYMMETRY
    # # center = IMP.algebra.Vector3D([174.619,175.188,177.687])
    # center = IMP.algebra.Vector3D([166.415,166.428,171.871])
    # rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    # transform = IMP.algebra.get_rotation_about_point(center, rot)
    #
    # for subunit in subunits:
    #     dof.constrain_symmetry(
    #         references=mols[subunit],
    #         clones=clones[subunit],
    #         transform=transform
    #     )

    m.update()

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

    # # SHUFFLE
    # shuffle_rbs = list()
    # shuffle_rbs.extend(rbs["MSIN1"][1:])
    # shuffle_rbs.extend(rbs["AKT1"])
    # IMP.pmi.tools.shuffle_configuration(
    #     objects=shuffle_rbs,
    #     max_translation=100,
    #     max_rotation=2*math.pi,
    #     verbose=False,
    #     cutoff=10,
    #     niterations=25
    # )
    # dof.optimize_flexible_beads(100)

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
    test_em_mtorc2(
        n_frames=int(sys.argv[1]),
        map=None if sys.argv[2] == "" else sys.argv[2]
    )