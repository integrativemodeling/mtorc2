import pandas as pd
from pathlib import Path
import IMP
import IMP.atom
import IMP.rmf
import RMF
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from rmf_to_pdb import rmf_to_pdb


'''
If a bead has resolution 10, the residue index of the bead (in a pdb or rmf file) is the average of the range of residues that the bead represents.
'''
def get_res_id_10(
        m,
        pid
):
    start, end = m.get_particle_name(pid).split("_")[0].split("-")
    start, end = int(start), int(end)
    bead_num = (end - start) / 2

    return start + math.ceil(bead_num)


def get_xl_vis_df(
        xl_df,
        m,
        h
):
    chain_ids = dict()
    # chain_ids["MTOR"] = "A"
    # chain_ids["RICTOR"] = "B"
    # chain_ids["MLST8"] = "C"
    # chain_ids["MSIN1"] = "D"
    # chain_ids["CRIM"] = "E"
    # chain_ids["RBD"] = "F"
    # chain_ids["MSIN1PH"] = "G"
    # chain_ids["AKT1PH"] = "H"
    # chain_ids["KINASE"] = "I"

    chain_ids["MTOR"] = "A"
    chain_ids["RICTOR"] = "B"
    chain_ids["MLST8"] = "C"
    chain_ids["MSIN1"] = "D"
    chain_ids["AKT1"] = "E"

    # fh = RMF.open_rmf_file_read_only(str(rmf_file))
    # m = IMP.Model()
    # h = IMP.rmf.create_hierarchies(fh, m)[0]

    xl_vis_dict = dict()
    xl_vis_dict["chain1"] = list()
    xl_vis_dict["res1"] = list()
    xl_vis_dict["chain2"] = list()
    xl_vis_dict["res2"] = list()
    xl_vis_dict["group"] = list()

    for i in range(len(xl_df)):
        prot1,res1,prot2,res2,xl_type,cutoff,xl_sat,min_dist,ambig,frame = xl_df.iloc[i]
        # print(prot1, res1, prot2, res2, type, cutoff, xl_sat, ambig)

        res_1_vis, res_2_vis = res1, res2
        sel1 = IMP.atom.Selection(h, molecule=prot1, residue_index=int(res1))
        pid1 = sel1.get_selected_particle_indexes()[0]
        if "bead" in m.get_particle_name(pid1):
            res_1_vis = get_res_id_10(m, pid1)

        sel2 = IMP.atom.Selection(h, molecule=prot2, residue_index=int(res2))
        pid2 = sel2.get_selected_particle_indexes()[0]
        if "bead" in m.get_particle_name(pid2):
            res_2_vis = get_res_id_10(m, pid2)

        chain_1 = chain_ids[prot1]
        if ambig:
            chain_2 = chr(ord(chain_ids[prot2]) + 5)
        else:
            chain_2 = chain_ids[prot2]

        xl_vis_dict["chain1"].append(chain_1)
        xl_vis_dict["res1"].append(res_1_vis)
        xl_vis_dict["chain2"].append(chain_2)
        xl_vis_dict["res2"].append(res_2_vis)

        xl_vis_dict["group"].append(xl_sat)

    xl_vis_df = pd.DataFrame(xl_vis_dict)
    return xl_vis_df


def draw_xls(
        xl_vis_df,
        pb_file
):
    groups_dict = dict()
    groups_dict[0], groups_dict[1], groups_dict[2] = "", "", ""
    for i in range(len(xl_vis_df)):
        chain_1, res_1, chain_2, res_2, group = xl_vis_df.iloc[i, 0:5]

        entry = "#1/{}:{} #1/{}:{}\n".format(chain_1, res_1, chain_2, res_2)
        print(entry)
        if chain_1 == chain_2 and res_1 == res_2:
            continue

        groups_dict[group] = groups_dict[group] + entry

    colors_dict = {0: "red", 1: "green", 2: "#FF9900"}
    f = open(pb_file, "w")
    for group in groups_dict.keys():
        color = colors_dict[group]
        f.write("; halfbond = false\n; color = {}\n; radius = 1\n; dashes = 0\n".format(color))
        f.write(groups_dict[group])
    f.close()


if __name__ == "__main__":
    cluster_dir = Path("/wynton/home/sali/mhancock/mtorc2/data/final_models/submission_2/126_0_2_3/cluster.0")

    rmf_file = Path(cluster_dir, "cluster_center_model.rmf3")
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(0))

    xl_file = Path(cluster_dir, "xls.csv")
    xl_df = pd.read_csv(xl_file)
    xl_df = xl_df.drop(columns=["Unnamed: 0"])

    exclude_xl_files = list()
    exclude_xl_files.append(Path(Path.home(), "mtorc2/data/xlms/csvs/a_site.csv"))
    exclude_xl_files.append(Path(Path.home(), "mtorc2/data/xlms/csvs/akt_tail.csv"))
    exclude_dfs = list()
    for exclude_xl_file in exclude_xl_files:
        exclude_dfs.append(pd.read_csv(exclude_xl_file))

    exclude_df = pd.concat(exclude_dfs)

    exclude_indices = list()
    for i in range(len(xl_df)):
        exclude = False
        prot1,res1,prot2,res2,xl_type,cutoff,sat,min_dist,ambig,frame = xl_df.iloc[i]

        print(prot1,res1,prot2,res2)
        # len(exclude_df[(exclude_df["prot1"] == prot1) & (exclude_df["res1"] == res1) & (exclude_df["prot2"] == prot2) & (exclude_df["res2"] == res2)]))

        # exclude = exclude_xl(prot1, prot2, res1, res2, exclude_df)
        if len(exclude_df[(exclude_df["prot1"] == prot1) & (exclude_df["res1"] == res1) & (exclude_df["prot2"] == prot2) & (exclude_df["res2"] == res2)]) > 0:
            exclude = True
        elif len(exclude_df[(exclude_df["prot2"] == prot1) & (exclude_df["res2"] == res1) & (exclude_df["prot1"] == prot2) & (exclude_df["res1"] == res2)]) > 0:
            exclude = True

        if exclude:
            exclude_indices.append(i)

    xl_df = xl_df.drop(exclude_indices)
    print(len(xl_df))

    xl_vis_df = get_xl_vis_df(
        xl_df=xl_df,
        m=m,
        h=h
    )

    pb_file = Path(cluster_dir, "xls_new.pb")
    draw_xls(
        xl_vis_df=xl_vis_df,
        pb_file=pb_file
    )

