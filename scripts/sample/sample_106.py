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

def sample_106(
        output_dir,
        n_frames,
        em_comps,
        map,
        res_per_comp,
):
    print("output_dir:      {}".format(output_dir))
    print("n_frames:        {}".format(n_frames))
    print("res_per_comp     {}".format(res_per_comp))
    print("em_comps:        {}".format(em_comps))
    print("map              {}".format(map))

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.domain.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    structures = dict()
    pdb_file = Path(Path.home(), "mtorc2/data/pdb/mTORC2.domain.pdb")

    if "CRIM" in em_comps.split(","):
        structures["CRIM"] = (pdb_file, "E", "purple", "N")
    if "KINASE" in em_comps.split(","):
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

                if res_per_comp == 1:
                    mol.add_representation(
                        residues=atomic_residues,
                        # density_prefix=density_prefix,
                        density_force_compute=True,
                        resolutions=[1],
                        setup_particles_as_densities=True,
                        color=color
                    )
                elif res_per_comp > 1:
                    mol.add_representation(
                        residues=atomic_residues,
                        density_residues_per_component=res_per_comp,
                        density_prefix=density_prefix,
                        density_force_compute=False,
                        resolutions=[1,10],
                        color=color
                    )
            else:
                mol.add_representation(
                    residues=atomic_residues,
                    resolutions=[1,10],
                    color=color
                )

        flex = list(mol.get_non_atomic_residues())
        print(len(flex))

        print("REPRESENTATION")
        print(len(mol.get_residues()))
        print(len(flex))
        print(len(atomic_residues))

        setup_particles_as_densities = False
        print(component, setup_particles_as_densities)

        mol.add_representation(
            flex,
            resolutions=[10],
            color=color,
            setup_particles_as_densities=setup_particles_as_densities
        )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    # RIGID BODIES
    for component in mols.keys():
        mol_rb = create_rigid_body(
            mol=mols[component],
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

    m.update()

    # EM
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
    sample_106(
        output_dir=sys.argv[1],
        n_frames=int(sys.argv[2]),
        em_comps=sys.argv[3],
        map=sys.argv[4],
        res_per_comp=int(sys.argv[5]),
    )