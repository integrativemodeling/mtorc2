import pandas as pd
from pathlib import Path
import IMP
import IMP.atom
import IMP.rmf
import RMF
import random
import numpy as np
import matplotlib.pyplot as plt
import IMP.pmi.tools
from operator import itemgetter
import math
import time
import multiprocessing


def get_ps(rmf_file, frame_num):
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(frame_num))

    s = IMP.atom.Selection(h, chain_ids=["A", "B", "C", "D", "E"], resolution=10)
    pids = s.get_selected_particle_indexes()

    return h, m, pids


def calculate_distances(
        pool_param
):
    rmf_file = pool_param["rmf_file"]
    frame_num = pool_param["frame_num"]
    ref_rmf = pool_param["ref_rmf"]

    print(rmf_file, frame_num)

    h, m, pids = get_ps(rmf_file, frame_num)
    h_0, m_0, pids_0 = get_ps(ref_rmf, 0)
    s = IMP.atom.Selection(h, chain_id="A")
    s_0 = IMP.atom.Selection(h_0, chain_id="A")

    t = IMP.atom.get_transformation_aligning_first_to_second(s, s_0)

    dist_matrix = np.ndarray(shape=(len(pids), len(pids_0)))
    for i in range(len(pids)):
        pid_0 = pids[i]
        xyz_0 = IMP.core.XYZR(m, pid_0).get_coordinates()
        xyz_0_t = t.get_transformed(xyz_0)

        for j in range(i,len(pids)):
            pid_1 = pids[j]
            xyz_1 = IMP.core.XYZR(m, pid_1).get_coordinates()
            xyz_1_t = t.get_transformed(xyz_1)

            d = (xyz_0_t - xyz_1_t).get_magnitude()

            dist_matrix[i,j] = d
            dist_matrix[j,i] = d

    return dist_matrix, rmf_file.stem, frame_num


if __name__ == "__main__":
    gsms_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/gsms/1000R")
    cluster = 0
    job_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/1")
    cluster_dir = Path(job_dir, "cluster.{}".format(cluster))

    ref_rmf = Path(cluster_dir, "cluster_center_model.rmf3")

    with open(Path(gsms_dir, "A.txt"), 'r') as fp:
        N_A = len(fp.readlines())
        print("N_A: {}".format(N_A))

    samp_clust_file = Path(job_dir, "cluster.{}.all.txt".format(cluster))
    samp_clust_df = pd.read_csv(samp_clust_file, header=None)

    print(samp_clust_file)
    print(len(samp_clust_df))

    # cluster.$.all.txt contains a list of numbers which correspond to frame ids from either the A.rmf3 file or the B.rmf3 file. If the frame id is > than the number of A.rmf3 frames, then the frame id is really a B.rmf3 frame (subtract N_A from the frame id).
    pool_params = list()
    for i in range(len(samp_clust_df)):
    # for i in range(100):
        params_dict = dict()
        frame_num = samp_clust_df.iloc[i,0]
        if frame_num < N_A:
            rmf_file = Path(gsms_dir, "A.rmf3")
        else:
            rmf_file = Path(gsms_dir, "B.rmf3")
            frame_num = frame_num - N_A

        params_dict["rmf_file"] = rmf_file
        params_dict["frame_num"] = frame_num
        params_dict["ref_rmf"] = ref_rmf

        pool_params.append(params_dict)
        # break

    print("# inputs: {}".format(len(pool_params)))
    print("CPUs: {}".format(multiprocessing.cpu_count()))
    # for pool_param in pool_params:
    #     calculate_distances(pool_param)
    #     break

    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    # Calculate the distances for all RMF frames relative to the median structure and save each to a .npy matrix.
    outputs = pool_obj.map(
        calculate_distances,
        pool_params
    )

    matrix_dir = Path("/wynton/group/sali/mhancock/mtorc2/tmp")
    for dist_matrix, rmf_file_name, rmf_id in outputs:
        matrix_file = Path(matrix_dir, "dist_{}{}.npy".format(rmf_file_name, rmf_id))
        np.save(matrix_file, dist_matrix)

    mat_files = [mat_file for mat_file in matrix_dir.glob("dist_*")]
    print(len(mat_files))
    dist_mat = np.load(mat_files[0])
    avg_mat = np.zeros(dist_mat.shape)

    # Read all the .npy matrices and compute the RMSF matrix.
    for mat_file in mat_files:
        print(mat_file)
        dist_mat = np.load(mat_file)
        avg_mat = avg_mat + dist_mat

        # break

    rmsf_file = Path(job_dir, "cluster.0.rmsf.npy")
    avg_file = Path(job_dir, "cluster.0.avg.npy")
    avg_mat = avg_mat / len(mat_files)
    np.save(avg_file, avg_mat)

    rmsf_mat = np.zeros(dist_mat.shape)

    for mat_file in mat_files:
        print(mat_file)
        dist_mat = np.load(mat_file)
        rmsf_mat = rmsf_mat + (dist_mat - avg_mat)**2

        # break

    rmsf_mat = rmsf_mat / len(mat_files)
    rmsf_mat = rmsf_mat**(1/2)

    np.save(rmsf_file, rmsf_mat)
