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


def calculate_distances(
    model
):
    rmf_file, frame_num = model

    print(rmf_file, frame_num)
    included_objects = list()
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(frame_num))
    for i in range(0,10,2):
        h_mol = h.get_children()[0].get_children()[i]
        included_objects.append(h_mol)

    hierarchies = IMP.pmi.tools.input_adaptor(
        included_objects,
        10,
        flatten=True
    )
    mdl = hierarchies[0].get_model()
    included_ps = [h.get_particle() for h in hierarchies]

    dist_matrix = np.ndarray(shape=(len(included_ps), len(included_ps)))

    for i in range(len(included_ps)):
        for j in range(len(included_ps)):
            pt_0 = included_ps[i]
            pt_1 = included_ps[j]

            xyz_0 = IMP.core.XYZR(pt_0)
            xyz_1 = IMP.core.XYZR(pt_1)
            d = IMP.algebra.get_distance(xyz_0.get_sphere(), xyz_1.get_sphere())

            dist_matrix[i,j] = d
            dist_matrix[j,i] = d

    return dist_matrix


if __name__ == "__main__":
    gsms_dir = Path("/wynton/group/sali/mhancock/mtorc2/exp_13/119/analysis/6/sampcon_2/gsms/1000R")
    cluster = 0
    sampcon_job_dir = Path(Path.home(), "mtorc2/manuscript/submission/models/119_6_2_2")
    cluster_dir = Path(sampcon_job_dir, "cluster.{}".format(cluster))

    with open(Path(gsms_dir, "A.txt"), 'r') as fp:
        N_A = len(fp.readlines())
        print("N_A: {}".format(N_A))

    samp_clust_file = Path(sampcon_job_dir, "cluster.{}.all.txt".format(cluster))
    samp_clust_df = pd.read_csv(samp_clust_file, header=None)

    print(samp_clust_file)
    print(len(samp_clust_df))

    models = list()
    for i in range(len(samp_clust_df)):
        frame_num = samp_clust_df.iloc[i,0]
        if frame_num < N_A:
            rmf_file = Path(gsms_dir, "A.rmf3")
        else:
            rmf_file = Path(gsms_dir, "B.rmf3")
            frame_num = frame_num - N_A

        models.append((rmf_file, frame_num))

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    dist_matrices = pool_obj.map(
        calculate_distances,
        models
    )

    dist_matrix_shape = (dist_matrices[0].shape[0], dist_matrices[0].shape[1], len(models))
    dist_matrix_merge = np.ndarray(shape=dist_matrix_shape)

    for i in range(len(dist_matrices)):
        dist_matrix_merge[:,:,i] = dist_matrices[i]

    matrix_file = Path(cluster_dir, "distances.npy")
    np.save(matrix_file, dist_matrix_merge)