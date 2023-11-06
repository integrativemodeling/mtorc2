import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import time
from pathlib import Path
import sys
sys.path.append("/wynton/home/sali/mhancock/mtorc2/src")
from representation import get_representation
from representation import create_rigid_bodies
from representation import setup_c2_symmetry
from restraints import get_ev_restraints
from restraints import get_xl_restraint
from restraints import get_crs


def dimer_flex(
        xls,
        n_steps
):
    output_dir = Path("./output_0")

    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    mols, clones = get_representation(
        m=m,
        s=s,
        dimer=True,
        flex=True
    )

    root_hier = s.build()
    dof = IMP.pmi.dof.DegreesOfFreedom(m)

    create_rigid_bodies(
        dof=dof,
        mols=mols,
        clones=clones
    )

    output_objects = []  # keep a list of functions that need to be reported
    crs = get_crs(
        mols=mols
    )
    for cr in crs:
        cr.add_to_model()
        output_objects.append(cr)

    ev_rs = get_ev_restraints(
        mols=mols,
        clones=clones
    )
    for ev_r in ev_rs:
        ev_r.add_to_model()
        output_objects.append(ev_r)

    setup_c2_symmetry(
        mols=mols,
        clones=clones,
        dof=dof
    )
    m.update()

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

    # Don't shuffle components
    t0 = time.time()
    rex = IMP.pmi.macros.ReplicaExchange0(
        m,
        root_hier=root_hier,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=str(output_dir),
        output_objects=output_objects,
        monte_carlo_steps=10,
        number_of_best_scoring_models=0,
        number_of_frames=n_steps,
        replica_exchange_swap=False
    )

    rex.execute_macro()
    return time.time() - t0