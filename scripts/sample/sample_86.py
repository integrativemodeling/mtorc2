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


def create_rigid_body(
        mol,
        component,
        dof
):
    rbs = list()
    if component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
        max_trans, max_rot = 0,0
    else:
        max_trans, max_rot = 5,5

    print(component, max_trans, max_rot)

    if component == "KINASE":
        rigid_parts_list = list()
        # The first entry is the N lobe and the second entry is the C lobe.
        for start, end in [("1", "89"), ("90", "335")]:
            rigid_parts_list.append(
                mol.residue_range(
                    a=start,
                    b=end
                )
            )

        non_rigid_parts_list = list()
        # The first entry is the short loop connecting the N and C lobes. The second entry is the tail region of the C lobe (which may or not be ommitted from the structure).
        for start, end in [("84", "89"), ("273", "335")]:
            non_rigid_parts_list.append(
                mol.residue_range(
                    a=start,
                    b=end
                )
            )
        # non_rigid_parts_list.append(None)

    else:
        rigid_parts_list = [mol.get_residues()]
        non_rigid_parts_list = [mol.get_non_atomic_residues()]

    for i in range(len(rigid_parts_list)):
        rigid_parts = rigid_parts_list[i]
        non_rigid_parts = non_rigid_parts_list[i]

        print("RIGID BODIES")
        print(len(rigid_parts))
        print(len(non_rigid_parts))

        rb_movers, rb = dof.create_rigid_body(
            rigid_parts=rigid_parts,
            nonrigid_parts=non_rigid_parts,
            max_trans=max_trans,
            max_rot=max_rot,
            nonrigid_max_trans=5
        )
        rbs.append(rb)

    return rbs

def run_86(
        output_dir,
        n_frames,
        xls,
        conn,
        map,
        em_comps,
        res_per_comp,
        removed_regions,
        akt1_ph
):
    print("output_dir:      {}".format(output_dir))
    print("xls:             {}".format(xls))
    print("n_frames:        {}".format(n_frames))
    print("conn:            {}".format(conn))
    print("em_comps:        {}".format(em_comps))
    print("res_per_comp     {}".format(res_per_comp))
    print("map              {}".format(map))
    print("removed_regions: {}".format(removed_regions))
    print("akt1_ph:         {}".format(akt1_ph))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.domain.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(Path.home(), "mtorc2/data/pdb/mTORC2.domain.pdb")
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

        if component == "KINASE":
            res_ranges = [("N",(1,83)),("C",(90,272))]
        else:
            res_ranges = [("",(1,len(mol.get_residues())))]
        for termini, res_range in res_ranges:
            atomic_residues = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=chain,
                soft_check=True,
                res_range=res_range
            )

            if component in em_comps:
                density_prefix = component+termini
                density_residues_per_component = res_per_comp
            else:
                density_prefix = None
                density_residues_per_component = None

            print(density_prefix, density_residues_per_component)
            mol.add_representation(
                residues=atomic_residues,
                density_residues_per_component=density_residues_per_component,
                density_prefix=density_prefix,
                density_force_compute=False,
                resolutions=[1,10],
                color=color
            )

        flex = list(mol.get_non_atomic_residues())
        print(len(flex))
        for region in removed_regions.split(","):
            region_component = region.split(":")[0]
            if component == region_component:
                region_range = region.split(":")[1].split("-")
                start, end = region_range[0], region_range[1]
                # Strings for pdb indexing.
                res_range = mol.residue_range(
                    a=start,
                    b=end
                )
                flex = list(set(flex) - set(res_range))

        print("REPRESENTATION")
        print(len(mol.get_residues()))
        print(len(flex))
        print(len(atomic_residues))

        setup_particles_as_densities = False
        if component+"FLEX" in em_comps:
            setup_particles_as_densities = True
            print(component, setup_particles_as_densities)

        if component == "RICTOR" and "RICTORFLEX" in em_comps:
            mol.add_representation(
                flex[:25],
                resolutions=[10],
                color=color,
                setup_particles_as_densities=True
            )

            mol.add_representation(
                flex[25:],
                resolutions=[10],
                color=color,
                setup_particles_as_densities=False
            )
        elif len(flex) > 0:
                mol.add_representation(
                    flex,
                    resolutions=[10],
                    color=color,
                    setup_particles_as_densities=setup_particles_as_densities
                )

        clone = mol.create_clone(
            chain_id=clone_chain
        )
        clones[component] = clone

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    for component in mols.keys():
        mol_rb = create_rigid_body(
            mol=mols[component],
            component=component,
            dof=dof
        )

        clone_rb = create_rigid_body(
            mol=clones[component],
            component=component,
            dof=dof
        )

    output_objects = []
    # CONNECTIVITY
    for component in mols.keys():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            objects=mols[component],
            resolution=10
        )
        cr.add_to_model()
        output_objects.append(cr)

    def get_res_id_n_and_c(
            component
    ):
        res_id_n = int(str(list(mols[component].get_residues())[0]).split("_")[3][1:])
        res_id_c = int(str(list(mols[component].get_residues())[-1]).split("_")[3][1:])
        return res_id_n, res_id_c

    connections = list()
    connections.append(("MSIN1", "CRIM"))
    connections.append(("CRIM", "RBD"))
    connections.append(("RBD", "MSIN1PH"))
    for connection in connections:
        prot1, prot2 = connection
        res1, res2 = get_res_id_n_and_c(prot1)[1], get_res_id_n_and_c(prot2)[0]
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
        print(component)
        dof.constrain_symmetry(
            references=mols[component],
            clones=clones[component],
            transform=transform
        )

    m.update()

    # XLS
    for xl in xls:
        print(xl)
        xl, xl_w = xl.split(":")
        xl_w = int(xl_w)
        if xl == "":
            continue
        # file, weight, type = xls[i]
        xl_r = get_xl_restraint(
            hier=root_hier,
            xl_file=xl+".csv",
            w=xl_w,
            type="DSS"
        )
        xl_r.add_to_model()
        output_objects.append(xl_r)

    # EM
    if map:
        sel = IMP.atom.Selection(
            hierarchy=root_hier,
            representation_type=IMP.atom.DENSITIES
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

    print("Shuffling")
    shuffle_objs = [mols["KINASE"], clones["KINASE"]]
    IMP.pmi.tools.shuffle_configuration(
        objects=shuffle_objs,
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
        score_moved=False
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)


if __name__ == "__main__":
    run_86(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
        xls=sys.argv[3].split(","),
        conn=int(sys.argv[4]),
        em_comps=sys.argv[5].split(","),
        res_per_comp=int(sys.argv[6]),
        map=None if sys.argv[7] == "" else sys.argv[7],
        removed_regions=sys.argv[8],
        akt1_ph=sys.argv[9]
    )