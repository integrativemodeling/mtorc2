import pandas as pd
from pathlib import Path
import IMP
import IMP.atom
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import sys

import xl_satisfaction


def get_xl_result_dfs(
        cluster_dir,
        gsms_dir,
        xl_df
):
    # Returns a list containing tuple-elements: (rmf_file, frame_num)
    gsms_files = xl_satisfaction.get_cluster_gsms(
        cluster_dir=cluster_dir,
        gsms_dir=gsms_dir
    )

    pool_inputs = list()
    for frame_num, rmf_file in gsms_files:
        pool_inputs.append((frame_num, rmf_file, xl_df.copy()))

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())

    # for pool_input in pool_inputs:
    #     score_results = [xl_satisfaction.pool_f(pool_input)]

    #     break

    score_results = pool_obj.imap(
        xl_satisfaction.pool_f,
        pool_inputs
    )

    xl_dfs = list()
    for frame, xl_result_df in score_results:
        print(frame, len(xl_result_df))
        xl_dfs.append(xl_result_df)

    pool_obj.close()

    return xl_dfs


if __name__ == "__main__":
    xl_files = list()
    xl_files.append((Path(Path.home(), "mtorc2/data/xlms/csvs/dss.csv"), "dss", 35))
    xl_files.append((Path(Path.home(), "mtorc2/data/xlms/csvs/edc.csv"), "edc", 16))

    sampcon_job_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/1")
    cluster_dir = Path(sampcon_job_dir, "cluster.0")

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
    xl_df.reset_index(inplace=True,drop=True)

    for i in range(len(xl_df)):
        prot1, res1, prot2, res2, type, cutoff = xl_df.iloc[i,0:6]
        print(prot1, res1, prot2, res2, type, cutoff)

    xl_result_dfs = get_xl_result_dfs(
        cluster_dir=cluster_dir,
        gsms_dir=gsms_dir,
        xl_df=xl_df
    )

    # Subtract 10 for a_site XLs.
    sats = list()
    for xl_result_df in xl_result_dfs:
        sats.append(xl_result_df["sat"].sum() / (len(xl_result_df)-10))

    print(np.max(sats), np.min(sats))