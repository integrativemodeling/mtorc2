from pathlib import Path
from restraints import get_xl_restraint
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.algebra
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time


def run_58(
        output_dir,
        xls,
        max_dr,
        em_weight,
        n_frames,
        n_pdbs=1
):
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

        if component in ["MSIN1", "CRIM", "RBD", "MSIN1PH", "AKT1PH", "KINASE"]:
            # density_prefix = Path(Path.home(), "mtorc2/data/em", component)
            density_prefix = component
            mol.add_representation(
                residues=atom,
                density_residues_per_component=5,
                density_prefix=density_prefix,
                # density_force_compute=True,
                resolutions=[1,10],
                color=color
            )

            particles_as_densities = False
            if component in ["MSIN1"]:
                particles_as_densities = True

            flex = mol.get_non_atomic_residues()
            mol.add_representation(
                flex,
                resolutions=[10],
                setup_particles_as_densities=particles_as_densities,
                color=color
            )

        else:
            mol.add_representation(
                residues=atom,
                resolutions=[1,10],
                color=color
            )

            if component in ["RICTOR"]:
                flex = list(mol.get_residues())

                flex_selection = flex[0:25]
                flex_selection.extend(flex[1529:1604])

                mol.add_representation(
                    flex_selection,
                    resolutions=[10],
                    color=color
                )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    for component in mols.keys():
        if component in ["MTOR", "RICTOR", "MLST8", "MSIN1"]:
            max_trans, max_rot = 0,0
        elif component in ["CRIM", "KINASE"]:
            max_trans, max_rot = 0,5
        else:
            max_trans, max_rot = 5,5

        dof.create_rigid_body(
            rigid_parts=mols[component],
            nonrigid_parts=mols[component].get_non_atomic_residues(),
            max_trans=max_trans,
            max_rot=max_rot,
            nonrigid_max_trans=5
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

    if max_dr:
        connections = list()
        connections.append(("MSIN1", 157, "CRIM", 1))
        connections.append(("CRIM", 122, "RBD", 1))
        connections.append(("RBD", 102, "MSIN1PH", 1))
        connections.append(("AKT1PH", 145, "KINASE", 1))
        for connection in connections:
            prot1, res1, prot2, res2 = connection
            dr = IMP.pmi.restraints.basic.DistanceRestraint(
                root_hier=s.get_hierarchy(),
                tuple_selection1=(res1, res1, prot1, 0),
                tuple_selection2=(res2, res2, prot2, 0),
                distancemin=0,
                distancemax=max_dr,
                resolution=10,
                weight=1
            )
            dr.add_to_model()
            output_objects.append(dr)

    # EXCLUDED VOLUME
    ev_objects = list()
    for component in mols.keys():
        ev_objects.append(mols[component])

    ev_r = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=ev_objects,
        resolution=10
    )
    ev_r.add_to_model()
    output_objects.append(ev_r)

    m.update()

    # XLS
    for i in range(len(xls)):
        file, weight, type = xls[i]
        xl_r = get_xl_restraint(
            hier=root_hier,
            xl_file=file,
            w=weight,
            type=type
        )
        xl_r.add_to_model()
        output_objects.append(xl_r)

    # EM
    if em_weight:
        target_gmm_file = Path(glob_data_dir, "em/mtorc2_dsfact2_ng100.txt")
        em_components = IMP.pmi.tools.get_densities(root_hier)
        gemt = IMP.pmi.restraints.em.GaussianEMRestraint(
            em_components,
            str(target_gmm_file),
            # scale_target_to_mass=True,
            slope=0.000001,
            weight=em_weight
        )
        gemt.add_to_model()
        output_objects.append(gemt)

    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=20,
        number_of_best_scoring_models=n_pdbs,
        number_of_frames=n_frames
        # test_mode=True
    )

    t0 = time.time()
    rex.execute_macro()
    print((time.time() - t0)/n_frames)