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
from restraints import get_xl_restraint


def run_71(
        output_dir,
        xls,
        n_frames,
        conn,
        em,
        em_comps,
        akt_dr,
        akt1_ph,
        shuffle,
        partial_flex,
        move_aware,
        map
):
    print("output_dir:   {}".format(output_dir))
    print("xls:          {}".format(xls))
    print("n_frames:     {}".format(n_frames))
    print("conn:         {}".format(conn))
    print("em:           {}".format(em))
    print("em_comps:     {}".format(em_comps))
    print("akt_ph:       {}".format(akt1_ph))
    print("shuffle:      {}".format(shuffle))
    print("partial_flex: {}".format(partial_flex))
    print("move_aware:   {}".format(move_aware))
    print("map           {}".format(map))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.domain.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(Path.home(), "mtorc2/data/pdb/1601720.no_bead.pdb")
    structures["MTOR"] = (pdb_file, "A", "orange", "J")
    structures["RICTOR"] = (pdb_file, "B", "blue", "K")
    structures["MLST8"] = (pdb_file, "C", "green", "L")
    structures["MSIN1"] = (pdb_file, "D", "red", "M")
    structures["CRIM"] = (pdb_file, "E", "purple", "N")
    structures["RBD"] = (pdb_file, "F", "brown", "O")
    structures["MSIN1PH"] = (pdb_file, "G", "tan", "P")
    if akt1_ph:
        structures["AKT1PH"] = (pdb_file, "H", "salmon", "Q")
    structures["KINASE"] = (pdb_file, "I", "black", "R")

    clones = dict()
    mols = dict()
    for component in structures.keys():
        print(component)
        pdb_file, chain, color, clone_chain = structures[component]

        mol = st.create_molecule(
            name=component,
            sequence=seqs[component],
            chain_id=chain
        )
        mols[component] = mol

        atom = mol.add_structure(
            pdb_fn=str(pdb_file),
            chain_id=chain,
            soft_check=True
        )

        flex = list(mol.get_non_atomic_residues())

        if component in em_comps:
            density_prefix = component
            mol.add_representation(
                residues=atom,
                density_residues_per_component=10,
                density_prefix=density_prefix,
                density_force_compute=False,
                resolutions=[1,10],
                color=color
            )

            mol.add_representation(
                flex,
                resolutions=[10],
                color=color
            )

        else:
            mol.add_representation(
                residues=atom,
                resolutions=[1,10],
                color=color
            )

            if not partial_flex or component not in ["MLST8", "RBD", "KINASE"]:
                mol.add_representation(
                    flex,
                    resolutions=[10],
                    color=color
                )

        clone = mol.create_clone(
            chain_id=clone_chain
        )
        clones[component] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    rbs = list()
    for component in mols.keys():
        if component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
            max_trans, max_rot = 0,0
        else:
            max_trans, max_rot = 5,5

        print(component, max_trans, max_rot)
        rb_movers, rb = dof.create_rigid_body(
            rigid_parts=mols[component],
            nonrigid_parts=mols[component].get_non_atomic_residues(),
            max_trans=max_trans,
            max_rot=max_rot,
            nonrigid_max_trans=5
        )

        if component not in ["KINASE"]:
            rbs.append(rb)

    for component in clones.keys():
        if component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
            max_trans, max_rot = 0,0
        else:
            max_trans, max_rot = 5,5

        rb_movers, rb = dof.create_rigid_body(
            rigid_parts=clones[component],
            nonrigid_parts=clones[component].get_non_atomic_residues(),
            max_trans=max_trans,
            max_rot=max_rot,
            nonrigid_max_trans=5
        )

        if component not in ["KINASE"]:
            rbs.append(rb)

    output_objects = []
    # CONNECTIVITY
    for component in mols.keys():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=10
        )
        cr.add_to_model()
        output_objects.append(cr)

    connections = list()
    connections.append(("MSIN1", 157, "CRIM", 1))
    connections.append(("CRIM", 122, "RBD", 1))

    if partial_flex:
        connections.append(("RBD", 74, "MSIN1PH", 1))
    else:
        connections.append(("RBD", 102, "MSIN1PH", 1))

    if akt1_ph:
        connections.append(("AKT1PH", 145, "KINASE", 1))

    for connection in connections:
        prot1, res1, prot2, res2 = connection
        dr = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=s.get_hierarchy(),
            tuple_selection1=(res1, res1, prot1, 0),
            tuple_selection2=(res2, res2, prot2, 0),
            distancemin=0,
            distancemax=conn,
            resolution=10,
            weight=1
        )
        dr.add_to_model()
        output_objects.append(dr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for component in mols.keys():
        ev_objects.append(mols[component])
    for component in clones.keys():
        ev_objects.append(clones[component])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    # SYMMETRY
    center = IMP.algebra.Vector3D([174.619,175.188,177.687])
    rot = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0], math.pi)
    transform = IMP.algebra.get_rotation_about_point(center, rot)

    for component in clones.keys():
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    m.update()

    # XLS
    for xl in xls:
        print(xl)
        if xl == "":
            continue
        # file, weight, type = xls[i]
        xl_r = get_xl_restraint(
            hier=root_hier,
            xl_file=xl+".csv",
            w=1,
            type="DSS"
        )
        xl_r.add_to_model()
        output_objects.append(xl_r)

    # EM
    if em:
        sel = IMP.atom.Selection(
            hierarchy=root_hier,
            representation_type=IMP.atom.DENSITIES
        )
        densities = sel.get_selected_particles()
        target_gmm_file = Path(glob_data_dir, "em/maps/{}_dsfact2_ng100.txt".format(map))
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

    if akt_dr:
        akt_xl_df = pd.read_csv(Path(Path.home(), "mtorc2/data/xlms/akt.csv"))
        akt_xl_df = akt_xl_df.drop(columns=["Unnamed: 0"])
        for i in range(len(akt_xl_df)):
            prot1, res1, prot2, res2 = tuple(akt_xl_df.iloc[i])
            print(prot1, res1, prot2, res2)
            d_1 = IMP.pmi.restraints.basic.DistanceRestraint(
                root_hier=s.get_hierarchy(),
                tuple_selection1=(res1, res1, prot1, 0),
                tuple_selection2=(res2, res2, prot2, 0),
                distancemin=0,
                distancemax=35,
                resolution=10,
                weight=1
            )
            d_1.add_to_model()
            output_objects.append(d_1)

    print(len(dof.get_rigid_bodies()))
    print(dof.get_rigid_bodies())
    print(len(rbs))
    print(rbs)
    if shuffle:
        shuffle_objs = [mols["KINASE"], clones["KINASE"]]
        if akt1_ph:
            shuffle_objs.append(mols["AKT1PH"])
            shuffle_objs.append(clones["AKT1PH"])

        print("Shuffling")
        IMP.pmi.tools.shuffle_configuration(
            objects=shuffle_objs,
            # objects=root_hier,
            # excluded_rigid_bodies=dof.get_rigid_bodies(),
            max_translation=300,
            max_rotation=2*math.pi,
            verbose=False,
            cutoff=10,
            niterations=100
        )
        dof.optimize_flexible_beads(100)

    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=1,
        number_of_frames=n_frames,
        replica_exchange_maximum_temperature=5,
        score_moved=bool(move_aware)
    )

    t0 = time.time()
    rex.execute_macro()
    print("move_aware: {}".format(move_aware))
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    run_71(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
        xls=sys.argv[3].split(","),
        conn=int(sys.argv[4]),
        em=int(sys.argv[5]),
        em_comps=sys.argv[6].split(","),
        akt_dr=int(sys.argv[7]),
        akt1_ph=int(sys.argv[8]),
        shuffle=int(sys.argv[9]),
        partial_flex=int(sys.argv[10]),
        move_aware=int(sys.argv[11]),
        map=sys.argv[12]
    )