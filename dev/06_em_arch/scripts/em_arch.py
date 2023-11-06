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


def em_arch(
        output_dir,
        n_frames,
        map,
        ev_rbs,
        res_per_comp,
        crim_end
):
    print("output_dir:      {}".format(output_dir))
    print("n_frames:        {}".format(n_frames))
    print("map              {}".format(map))
    print("ev_rbs           {}".format(ev_rbs))
    print("res_per_comp:    {}".format(res_per_comp))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    subunits = ["MSIN1", "AKT1"]
    structures = dict()
    pdb_dir = Path(Path.home(), "mtorc2/data/pdb")
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_dir, "mtorc2_238.pdb"), "G", "red", 1, 157, 1, 157, 0, None))
    # structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, crim_end+116, -116, "CRIM"))
    structures["MSIN1"].append((Path(pdb_dir, "2rvk.pdb"), "A", "red", 158, 279, 274, crim_end+116, -116, "CRIM"))
    structures["AKT1"] = list()
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 145, 231, 145, 231, 0, "KINASEN"))
    structures["AKT1"].append((Path(pdb_dir, "3o96.pdb"), "A", "purple", 232, 480, 232, 429, 0, "KINASEC"))

    mols = dict()
    chains = dict()
    chains["MSIN1"] = "D"
    chains["AKT1"] = "E"

    # for subunit in subunits.keys():
    for subunit in subunits:
        print(subunit)
        # pdb_file, chain, color, clone_chain = structures[component]

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chains[subunit]
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

            if prefix:
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
            print(atom_res)

        if subunit == "MSIN1":
            flex_parts = mol.residue_range(
                a=str(133),
                b=str(157)
            )
        # elif subunit == "AKT1":
        #     flex_parts = mol.residue_range(
        #         a=str(229),
        #         b=str(234)
        #     )

            print(flex_parts)
            mol.add_representation(
                flex_parts,
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )

        print(subunit, len(mol.get_residues()))
        print(subunit, len(mol.get_atomic_residues()))
        print(subunit, len(mol.get_non_atomic_residues()))

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rb_ranges, rbs = dict(), dict()
    rb_ranges["MSIN1"] = [(1, 157, 0), (158, 279, 5)]
    rb_ranges["AKT1"] = [(145, 231, 5), (232, 480, 5)]

    for subunit in subunits:
        rbs[subunit] = list()
        mol = mols[subunit]
        for rb_start, rb_end, rb_trans in rb_ranges[subunit]:
            for copy in [mol]:
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

    output_objects = []
    # # CONNECTIVITY
    cr_main_comps = ["MSIN1", "AKT1"]

    for component in cr_main_comps:
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=1
        )
        cr.add_to_model()
        output_objects.append(cr)


#
    # cr_range = mols["MSIN1"].residue_range(
    #     a=str(1),
    #     b=str(279)
    # )
    # cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
    #     objects=cr_range,
    #     resolution=1
    # )
    # cr.add_to_model()
    # output_objects.append(cr)

    # EXCLUDED VOLUME
    # for component in mols.keys():
    rb_ranges["MSIN1"] = [(1, 157, 0), (158, 279, 5)]
    rb_ranges["AKT1"] = [(145, 234, 5), (235, 480, 5)]

    ev_objects = list()
    if "KinaseC" in ev_rbs:
        ev_objects.append(
            mols["AKT1"].residue_range(
                a=str(145),
                b=str(234)
            )
        )

    if "KinaseN" in ev_rbs:
        ev_objects.append(
            mols["AKT1"].residue_range(
                a=str(235),
                b=str(480)
            )
        )

    if "Crim" in ev_rbs:
        ev_objects.append(
            mols["MSIN1"].residue_range(
                a=str(158),
                b=str(279)
            )
        )
    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

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
        target_is_rigid_body=False
    )
    gem.add_to_model()
    gem.set_label("BayesianEM")
    output_objects.append(gem)

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
    em_arch(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
        map=None if sys.argv[3] == "" else sys.argv[3],
        ev_rbs=sys.argv[4].split(","),
        res_per_comp=int(sys.argv[5]),
        crim_end=int(sys.argv[6])
    )