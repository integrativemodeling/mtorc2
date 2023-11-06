import pandas as pd
from pathlib import Path
import IMP
import IMP.rmf
import RMF
import pandas as pd
import IMP.atom
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
import sys
sys.path.append("/home/matthew/mtorc2/single_traj/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/single_traj/src")
from analysis import *
import math


def xl_satisfaction(
        exp_num,
        traj_num,
        job_num,
        xls,
        symm,
):
    pdb_file = Path(Path.home(), "mtorc2/single_traj/exp_{}/{}/sample/output.{}/output_0/pdbs/model.0.pdb".format(exp_num, traj_num, job_num))
    print(pdb_file)
    xl_dir = Path(Path.home(), "mtorc2/data/xlms")
    xl_df = pd.DataFrame()
    for xl in xls:
        xl_file = Path(xl_dir, "{}.csv".format(xl))
        xl_df = pd.concat([xl_df, pd.read_csv(xl_file)])

    if len(xl_df) > 0:
        xl_df = xl_df.drop(columns=["Unnamed: 0"])

    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m)
    pids = list(m.get_particle_indexes())

    chain_ids, inv_chain_ids = dict(), dict()
    chain_ids["MTOR_0"] = "A"
    chain_ids["RICTOR_0"] = "B"
    chain_ids["MLST8_0"] = "C"
    chain_ids["MSIN1_0"] = "D"
    chain_ids["CRIM_0"] = "E"
    chain_ids["RBD_0"] = "F"
    chain_ids["MSIN1PH_0"] = "G"
    chain_ids["AKT1PH_0"] = "H"
    chain_ids["KINASE_0"] = "I"

    chain_ids["MTOR_1"] = "J"
    chain_ids["RICTOR_1"] = "K"
    chain_ids["MLST8_1"] = "L"
    chain_ids["MSIN1_1"] = "M"
    chain_ids["CRIM_1"] = "N"
    chain_ids["RBD_1"] = "O"
    chain_ids["MSIN1PH_1"] = "P"
    chain_ids["AKT1PH_1"] = "Q"
    chain_ids["KINASE_1"] = "R"

    for chain_id in chain_ids.keys():
        inv_chain_ids[chain_ids[chain_id]] = chain_id

    # Construct the dictionary of all atomic positions.
    atoms_dict = {"prot": [], "copy": [], "res": [], "xyz": [], "bead": []}
    for pid in pids:
        p_name = m.get_particle_name(pid)
        if IMP.atom.Residue.get_is_setup(m, pid):
            res = IMP.atom.Residue(m, pid)
            chain_id = IMP.atom.Chain(res.get_parent()).get_id()
            prot, copy = inv_chain_ids[chain_id].split("_")

            atoms_dict["prot"].append(prot)
            atoms_dict["copy"].append(int(copy))
            atoms_dict["res"].append(int(res.get_index()))
            atoms_dict["xyz"].append(IMP.core.XYZ(res.get_child(0)).get_coordinates())
            if str(res.get_residue_type()).split("\"")[1] == "BEA":
                atoms_dict["bead"].append(1)
            else:
                atoms_dict["bead"].append(0)

    atoms_df = pd.DataFrame(atoms_dict)

    chain_ids_inv = dict()
    for key in chain_ids.keys():
        chain_ids_inv[chain_ids[key]] = key

    contains = list()
    copies1 = list()
    copies2 = list()
    xl_sats = list()
    # Iterate through all xls and check if the xl is contained and satisfied.
    for i in range(len(xl_df)):
        prot1, res1, prot2, res2 = xl_df.iloc[i,0:4]

        coord1_0 = get_coords(atoms_df, prot1, res1, 0)
        coord1_1 = get_coords(atoms_df, prot1, res1, 1)
        coord2_0 = get_coords(atoms_df, prot2, res2, 0)
        coord2_1 = get_coords(atoms_df, prot2, res2, 1)

        contain = False
        xl_sat = False
        l2norm_0, l2norm_1 = 0,0
        copy1, copy2 = 0,0
        if coord1_0 and coord2_0:
            contain = True

            l2norm_0 = np.linalg.norm(np.array(coord1_0)-np.array(coord2_0))

            if symm:
                l2norm_1 = np.linalg.norm(np.array(coord1_0)-np.array(coord2_1))
            else:
                l2norm_1 = l2norm_0 + 1

            if (l2norm_0 < 35) or (l2norm_1 < 35):
                xl_sat = True

        if l2norm_1 < l2norm_0:
            copy2 = 1

        copies1.append(copy1)
        copies2.append(copy2)
        xl_sats.append(xl_sat)
        contains.append(contain)

        print(i, prot1, res1, prot2, res2, contain, xl_sat)

    xl_df["contains"] = contains
    xl_df["xl_sat"] = xl_sats
    xl_df["copy_1"] = copies1
    xl_df["copy_2"] = copies2

    sat_xls = str()
    unsat_xls = str()
    # Iterate through the satisfaction df and construct the ChimeraX psuedobond.
    for i in range(len(xl_df)):
        if xl_df.iloc[i, xl_df.columns.get_loc("contains")]:
            prot1 = xl_df.iloc[i, xl_df.columns.get_loc("prot1")]
            prot2 = xl_df.iloc[i, xl_df.columns.get_loc("prot2")]

            res1 = xl_df.iloc[i, xl_df.columns.get_loc("res1")]
            res2 = xl_df.iloc[i, xl_df.columns.get_loc("res2")]

            copy1 = xl_df.iloc[i, xl_df.columns.get_loc("copy_1")]
            copy2 = xl_df.iloc[i, xl_df.columns.get_loc("copy_2")]

            chain1 = chain_ids[xl_df.iloc[i, xl_df.columns.get_loc("prot1")]+"_"+str(copy1)]
            chain2 = chain_ids[xl_df.iloc[i, xl_df.columns.get_loc("prot2")]+"_"+str(copy2)]

            if is_res_10(atoms_df, prot1, res1, copy1):
                res1 = get_closest_bead(atoms_df, prot1, res1, copy1)

            if is_res_10(atoms_df, prot2, res2, copy2):
                res2 = get_closest_bead(atoms_df, prot2, res2, copy2)

            if prot1 == prot2 and res1 == res2 and copy1 == copy2:
                continue

            if xl_df.iloc[i, xl_df.columns.get_loc("xl_sat")]:
                sat_xls = sat_xls + "#1/{}:{}@CA #1/{}:{}@CA\n".format(chain1, res1, chain2, res2)
            else:
                unsat_xls = unsat_xls + "#1/{}:{}@CA #1/{}:{}@CA\n".format(chain1, res1, chain2, res2)

    pb_file = Path(Path.home(), "mtorc2/single_traj/exp_{}/{}/analysis/{}.pb".format(exp_num, traj_num, job_num))
    f = open(pb_file, "w")
    f.write("; halfbond = false\n; color = green\n; radius = 1\n; dashes = 0\n")
    f.write(sat_xls)
    f.write("; halfbond = false\n; color = red\n; radius = 1\n; dashes = 0\n")
    f.write(unsat_xls)
    f.close()

    # Create plots.
    n_contains = dict()
    n_satisfied = dict()
    n_xls = dict()
    prots = ["MTOR", "RICTOR", "MSIN1", "MLST8", "AKT1"]
    for prot in prots:
        for prot2 in prots:
            entry = prot + "|" + prot2
            if prot2 + "|" + prot not in n_contains.keys():
                n_contains[entry] = 0
                n_satisfied[entry] = 0
                n_xls[entry] = 0

    for i in range(len(xl_df)):
        prot1 = xl_df.iloc[i, xl_df.columns.get_loc('prot1')]
        prot2 = xl_df.iloc[i, xl_df.columns.get_loc('prot2')]

        if prot1 == "CRIM" or prot1 == "RBD" or prot1 == "MSIN1PH":
            prot1 = "MSIN1"
        if prot2 == "CRIM" or prot2 == "RBD" or prot2 == "MSIN1PH":
            prot2 = "MSIN1"

        if prot1 == "KINASE" or prot1 == "AKT1PH":
            prot1 = "AKT1"
        if prot2 == "KINASE" or prot2 == "AKT1PH":
            prot2 = "AKT1"

        entry = prot1 + "|" + prot2
        if entry not in n_contains.keys():
            entry = prot2 + "|" + prot1

        n_xls[entry] = n_xls[entry] + 1

        if xl_df.iloc[i, xl_df.columns.get_loc('contains')]:
            n_contains[entry] = n_contains[entry] + 1

        if xl_df.iloc[i, xl_df.columns.get_loc('xl_sat')]:
            n_satisfied[entry] = n_satisfied[entry] + 1

    fig, ax = plt.subplots(figsize=(25, 18), dpi=80)
    labels = list(n_contains.keys())
    contains_list = [n_contains[key] for key in n_contains.keys()]
    satisfied_list = [n_satisfied[key] for key in n_contains.keys()]
    n_xls_list = [n_xls[key] for key in n_contains.keys()]

    x = np.arange(len(labels))  # the label locations
    width = .75  # the width of the bars

    rects3 = ax.bar(x, n_xls_list, width/3, label='Total')
    rects1 = ax.bar(x + width/3, contains_list, width/3, label='Contains')
    rects2 = ax.bar(x + 2*width/3, satisfied_list, width/3, label='Satisfied')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Counts')
    ax.set_title('XL Counts')
    ax.set_xticks(x, labels)
    ax.set_xticklabels(labels, rotation=45)
    ax.legend()

    ax.bar_label(rects1, padding=1)
    ax.bar_label(rects2, padding=1)
    ax.bar_label(rects3, padding=1)

    fig.tight_layout()
    fig_file = Path(Path.home(), "mtorc2/single_traj/exp_{}/{}/analysis/{}.png".format(exp_num, traj_num, job_num))
    plt.savefig(fig_file)


if __name__ == "__main__":
    exp_infos = list()
    # exp_info.append(("exp_6/45", ["45.csv", "akt.csv"]))
    # exp_info.append(("exp_6/46", ["45.csv", "akt.csv"]))
    # exp_info.append(("exp_6/47", ["45.A.csv", "akt.csv"]))
    # exp_info.append(("exp_6/48", ["45.B.csv", "akt.csv"]))
    # exp_info.append(("exp_6/49", ["45.csv", "akt.csv"]))
    # exp_info.append(("exp_6/50", ["45.A.csv", "akt.csv"]))
    # exp_info.append(("exp_6/51", ["45.B.csv", "akt.csv"]))
    # exp_info.append(("exp_6/52", ["45.A.csv", "akt.csv"]))
    # exp_info.append(("exp_6/53", ["45.A.csv", "akt.csv"]))
    # exp_info.append(("exp_6/54", ["44.csv", "akt.csv"]))
    # exp_info.append(("exp_6/55", ["44.csv", "akt.csv"]))
    # exp_info.append(("exp_6/56", ["44.csv", "akt.csv"]))
    # exp_info.append(("exp_7/57", ["57.csv", "A.csv"], True))
    # exp_info.append(("exp_7/58", ["57.csv", "A.csv", "akt.csv"], False))
    # exp_info.append(("exp_7/62", ["57.csv", "A.csv", "akt.csv"], False))
    # exp_info.append((8, 65, 13385, [], False))
    exp_infos.append((8, 66, 81257, ["57", "A", "akt"], False))
    exp_infos.append((8, 67, 81258, ["57", "B", "akt"], False))
    exp_infos.append((8, 68, 81259, ["57", "A", "B", "akt"], False))
    exp_infos.append((8, 69, 81261, ["57", "A", "B", "akt"], False))
    exp_infos.append((8, 71, 81269, ["57", "A", "akt"], False))


for info in exp_infos:
        print(info)
        xl_satisfaction(
            exp_num=info[0],
            traj_num=info[1],
            job_num=info[2],
            xls=info[3],
            symm=info[4]
        )


