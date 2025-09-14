
import pandas as pd

import IMP
import IMP.atom


def get_atoms_df(m, pids, inv_chain_ids):
    atoms_dict = {"prot": [], "copy": [], "res": [], "xyz": [], "bead": []}
    for pid in pids:
        p_name = m.get_particle_name(pid)
        if IMP.atom.Residue.get_is_setup(m, pid):
            res = IMP.atom.Residue(m, pid)
            chain_id = IMP.atom.Chain(res.get_parent()).get_id()
            prot, copy = inv_chain_ids[chain_id].split("_")
            print(chain_id, prot, copy)

            atoms_dict["prot"].append(prot)
            atoms_dict["copy"].append(int(copy))
            atoms_dict["res"].append(int(res.get_index()))
            atoms_dict["xyz"].append(IMP.core.XYZ(res.get_child(0)).get_coordinates())
            if str(res.get_residue_type()).split("\"")[1] == "BEA":
                atoms_dict["bead"].append(1)
            else:
                atoms_dict["bead"].append(0)

    atoms_df = pd.DataFrame(atoms_dict)
    atoms_df.head(10)