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


def get_distances(
        m, h,
        prot_1,
        res_id_1,
        prot_2,
        res_id_2
):
    # print(prot_1, res_id_1, prot_2, res_id_2)
    # print(type(prot_1), type(res_id_1), type(prot_2), type(res_id_2))

    try:
        sel1 = IMP.atom.Selection(
            h,
            molecule=prot_1,
            residue_index=res_id_1,
            copy_index=0,
        )

        pid1 = sel1.get_selected_particle_indexes()[0]
        xyz1 = IMP.core.XYZ(m, pid1).get_coordinates()

        sel2_1 = IMP.atom.Selection(
            h,
            molecule=prot_2,
            residue_index=res_id_2,
            copy_index=0
        )

        pid2_1 = sel2_1.get_selected_particle_indexes()[0]
        xyz2_1 = IMP.core.XYZ(m, pid2_1).get_coordinates()

        sel2_2 = IMP.atom.Selection(
            h,
            molecule=prot_2,
            residue_index=res_id_2,
            copy_index=1
        )
        pid2_2 = sel2_2.get_selected_particle_indexes()[0]
        xyz2_2 = IMP.core.XYZ(m, pid2_2).get_coordinates()

        dist_1 = np.linalg.norm(np.array(xyz1)-np.array(xyz2_1))
        dist_2 = np.linalg.norm(np.array(xyz1)-np.array(xyz2_2))
    except IndexError:
        print("{} {}, {} {} FAIL".format(prot_1, res_id_1, prot_2, res_id_2))
        return 0, 0

    return dist_1, dist_2


def check_satisfaction(
        xl_df,
        m, h
):
    for i in range(len(xl_df)):
        prot1, res1, prot2, res2, type, cutoff = xl_df.iloc[i,0:6]

        # print(prot1,res1,prot2,res2)

        dist_1, dist_2 = get_distances(
            m, h,
            prot_1=prot1,
            res_id_1=res1,
            prot_2=prot2,
            res_id_2=res2
        )

        min_dist = dist_1
        if dist_2 < dist_1:
            min_dist = dist_2

        xl_df.iloc[i, xl_df.columns.get_loc("min_dist")] = min_dist
        if min_dist < cutoff:
            xl_df.iloc[i, xl_df.columns.get_loc("sat")] = 1

    return xl_df


def pool_f(
        pool_input
):
    rmf_file, rmf_frame, xl_df = pool_input
    print(rmf_file, rmf_frame)

    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(rmf_frame))

    xl_df = check_satisfaction(
        xl_df=xl_df,
        m=m, h=h
    )

    del fh

    return rmf_file.stem + str(rmf_frame), xl_df


def check_ambiguity(
        rmf_file,
        xl_df
):
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(0))

    for i in range(len(xl_df)):
        prot1, res1, prot2, res2 = xl_df.iloc[i,0:4]

        dist_1, dist_2 = get_distances(
            m, h,
            prot_1=prot1,
            res_id_1=res1,
            prot_2=prot2,
            res_id_2=res2
        )

        if dist_2 < dist_1:
            xl_df.iloc[i, xl_df.columns.get_loc("ambig")] = 1

    return xl_df


def get_cluster_gsms(
        cluster_dir,
        gsms_dir
):
    # Script needs to read in the number of structures contained in the A.rmf3 file because the structures in sampcon cluster A are indexed between 0:N_A-1 while structures in cluster B are indexed between N_A:N_A+N_B-1.
    with open(Path(gsms_dir, "A.txt"), 'r') as fp:
        N_A = len(fp.readlines())
        print("N_A: {}".format(N_A))

    # print(len(xl_df))
    # print(xl_df.head())

    samp_clust_file = Path(str(cluster_dir)+".all.txt")
    samp_clust_df = pd.read_csv(samp_clust_file, header=None)
    # samp_clust_df = samp_clust_df.iloc[:1]

    print(samp_clust_file)
    print(len(samp_clust_df))

    structure_files = list()
    for i in range(len(samp_clust_df)):
        frame_num = samp_clust_df.iloc[i,0]
        if frame_num < N_A:
            rmf_file = Path(gsms_dir, "A.rmf3")
        else:
            rmf_file = Path(gsms_dir, "B.rmf3")
            frame_num = frame_num - N_A

        structure_files.append((rmf_file, frame_num))

    return structure_files


def get_all_gsms(
        gsms_dir
):
    a_scores_file = Path(gsms_dir, "A.txt")
    b_scores_file = Path(gsms_dir, "B.txt")

    with open(a_scores_file, "r") as f:
        n_a = len(f.readlines())
    with open(b_scores_file, "r") as f:
        n_b = len(f.readlines())

    gsms_files = list()
    for i in range(n_a):
        rmf_file = Path(gsms_dir, "A.rmf3")
        gsms_files.append((rmf_file, i))
    for i in range(n_b):
        rmf_file = Path(gsms_dir, "B.rmf3")
        gsms_files.append((rmf_file, i))

    return gsms_files


def xl_satisfaction_cluster(
        sampcon_job_dir,
        cluster,
        gsms_dir,
        xl_df
):
    # Returns a list containing tuple-elements: (rmf_file, frame_num)
    if cluster == "all":
        gsms_files = get_all_gsms(
            gsms_dir=gsms_dir
        )
    else:
        gsms_files = get_cluster_gsms(
            cluster_dir=Path(sampcon_job_dir, "cluster.{}".format(cluster)),
            gsms_dir=gsms_dir
        )

    pool_inputs = list()
    for frame_num, rmf_file in gsms_files:
        pool_inputs.append((frame_num, rmf_file, xl_df.copy()))

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())

    # for pool_input in pool_inputs:
    #     score_results = [pool_f(pool_input)]
    #
    #     break

    score_results = pool_obj.map(
        pool_f,
        pool_inputs
    )

    for frame, xl_result_df in score_results:
        # unsat_xls = list(xl_df[xl_df["sat"] == 0].index)
        # for i in unsat_xls:
        #     if xl_result_df.iloc[i, xl_df.columns.get_loc("sat")]:
        #         xl_df.iloc[i, xl_df.columns.get_loc("sat")] = 1
        #         xl_df.iloc[i, xl_df.columns.get_loc("frame")] = rmf_frame
        #         # print()

        for i in range(len(xl_df)):
            min_dist = xl_result_df.iloc[i, xl_result_df.columns.get_loc("min_dist")]
            if min_dist < xl_df.iloc[i, xl_df.columns.get_loc("min_dist")]:
                xl_df.iloc[i, xl_df.columns.get_loc("min_dist")] = min_dist
                xl_df.iloc[i, xl_df.columns.get_loc("frame")] = frame

            if not xl_df.iloc[i, xl_df.columns.get_loc("sat")] and xl_result_df.iloc[i, xl_result_df.columns.get_loc("sat")]:
                xl_df.iloc[i, xl_df.columns.get_loc("sat")] = 1

    pool_obj.close()

    return xl_df


if __name__ == "__main__":
    # exp_dir = "exp_14/126"
    # analysis_id = "0"
    # hdbscan_cluster = "2"
    # sampcon_id = "3"
    # gsms_folder = "5000R"
    # clusters = [0]

    # home_dir = Path("/wynton/group/sali/mhancock")
    # sampcon_dir = Path(home_dir, "mtorc2/samples", exp_dir, "analysis", analysis_id, "sampcon_{}".format(hdbscan_cluster))
    # if sampcon_id:
    #     sampcon_job_dir = Path(sampcon_dir, sampcon_id)
    # else:
    #     sampcon_job_dir = None
    # gsms_dir = Path(sampcon_dir, "gsms", gsms_folder)

    xl_files = list()
    xl_files.append((Path(Path.home(), "mtorc2/data/xlms/csvs/dss.csv"), "dss", 35))
    xl_files.append((Path(Path.home(), "mtorc2/data/xlms/csvs/edc.csv"), "edc", 16))
    xl_files.append((Path(Path.home(), "mtorc2/data/xlms/csvs/edc_intra.csv"), "edc", 16))

    sampcon_job_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/1")
    cluster = 0
    cluster_dir = Path(sampcon_job_dir, "cluster.{}".format(cluster))

    gsms_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/gsms/1000R")

    frames = list()
    for file, type, cutoff in xl_files:
        tmp_df = pd.read_csv(file)
        tmp_df["type"] = type
        tmp_df["cutoff"] = cutoff
        frames.append(tmp_df)

    xl_df = pd.concat(frames)
    xl_df["sat"] = 0
    xl_df["min_dist"] = 1000
    xl_df["ambig"] = 0
    xl_df["frame"] = -1
    xl_df = xl_df.drop(columns=["Unnamed: 0"])

    rows = list()
    for i in range(len(xl_df)):
        prot1, res1, prot2, res2 = xl_df.iloc[i,0:4]
        rows.append(i)

    xl_df = xl_df.iloc[rows]
    xl_df.reset_index(
        inplace=True,
        drop=True
    )

    for i in range(len(xl_df)):
        prot1, res1, prot2, res2, type, cutoff = xl_df.iloc[i,0:6]
        # print(prot1, res1, prot2, res2, type, cutoff)

    xl_df = xl_satisfaction_cluster(
        sampcon_job_dir=sampcon_job_dir,
        cluster=cluster,
        gsms_dir=gsms_dir,
        xl_df=xl_df
    )

    if cluster_dir:
        rmf_file = Path(cluster_dir, "cluster_center_model.rmf3")
        xl_file = Path(cluster_dir, "xls_new.csv")

        # If testing against a cluster, need to check the centroid structure for ambiguous xl distances for visualization.
        xl_df = check_ambiguity(
            rmf_file=rmf_file,
            xl_df=xl_df
        )
    else:
        xl_file = Path(gsms_dir, "xls.csv")

    xl_df.to_csv(xl_file)