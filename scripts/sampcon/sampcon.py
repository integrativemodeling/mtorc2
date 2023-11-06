import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories
import shutil


def get_best_cluster(
        analysis_dir
):
    cluster_file = Path(analysis_dir, "traj/summary_hdbscan_clustering.dat")
    cluster_df = pd.read_csv(cluster_file)
    best_cluster = int(
        cluster_df[cluster_df.Total_Score == cluster_df.Total_Score.min()][
            "cluster"])

    return best_cluster


def save_n_rand_entries(
        orig_file,
        new_file,
        N
):
    orig_df = pd.read_csv(orig_file)
    if N > len(orig_df) or N == -1:
        new_df = orig_df
    else:
        random_entries = random.choices(list(orig_df.index), k=N)
        new_df = orig_df.iloc[random_entries].copy()

    new_df = new_df.rename(columns={'Unnamed: 0': ''})
    new_df.to_csv(new_file, index=False)


def save_n_best_entries(
        orig_file,
        new_file,
        N
):
    orig_df = pd.read_csv(orig_file)

    if N > len(orig_df) or N == -1:
        new_df = orig_df
    else:
        new_df = orig_df.nsmallest(
            n=N,
            columns="Total_Score"

        )

    new_df = new_df.rename(columns={'Unnamed: 0': ''})
    new_df.to_csv(new_file, index=False)
    print(new_df.head())


def extract_models(
        output_dir,
        analysis_dir,
        filter,
        gsms_dir,
        N,
        cluster
):
    print("**Running extract_models***")
    t0 = time.time()

    out_dirs = [str(out_dir) for out_dir in output_dir.glob("*") if
                out_dir.is_dir()]

    # gsms_dir = Path(sampcon_dir, "gsms")
    hdbscan_A_file = Path(analysis_dir, "traj/A_cluster{}.csv".format(cluster))
    hdbscan_B_file = Path(analysis_dir, "traj/B_cluster{}.csv".format(cluster))
    gsms_A_file = Path(gsms_dir, "A_cluster{}_{}.csv".format(cluster, N))
    gsms_B_file = Path(gsms_dir, "B_cluster{}_{}.csv".format(cluster, N))

    if filter == "R":
        save_n_entries = save_n_rand_entries
    elif filter == "B":
        save_n_entries = save_n_best_entries
    else:
        raise RuntimeError("INVALID FILTER PARAM")

    save_n_entries(
        orig_file=hdbscan_A_file,
        new_file=gsms_A_file,
        N=N
    )

    save_n_entries(
        orig_file=hdbscan_B_file,
        new_file=gsms_B_file,
        N=N
    )

    # Load module
    AT = AnalysisTrajectories(
        out_dirs,
        dir_name='output_',
        analysis_dir=str(Path(analysis_dir, "traj")),
        nproc=50
    )

    HA = AT.get_models_to_extract(str(gsms_A_file))
    HB = AT.get_models_to_extract(str(gsms_B_file))

    rmf_A_file = Path(gsms_dir, "A.rmf3")
    rmf_B_file = Path(gsms_dir, "B.rmf3")
    scores_A_file = Path(gsms_dir, "A.txt")
    scores_B_file = Path(gsms_dir, "B.txt")

    AT.extract_models_to_single_rmf(
        HA,
        str(rmf_A_file),
        str(output_dir),
        str(scores_A_file)
    )

    AT.extract_models_to_single_rmf(
        HB,
        str(rmf_B_file),
        str(output_dir),
        str(scores_B_file)
    )

    print("Finished in {}s".format(time.time()-t0))


def sampcon(
        analysis_dir,
        gsms_dir,
        mode,
        thresh,
        range
):
    print("**Running sampcon***")
    t0 = time.time()

    rmf_A_file = Path(gsms_dir, "A.rmf3")
    rmf_B_file = Path(gsms_dir, "B.rmf3")
    scores_A_file = Path(gsms_dir, "A.txt")
    scores_B_file = Path(gsms_dir, "B.txt")

    custom_ranges_dir = Path(Path.home(), "mtorc2/data/custom_ranges")
    selection_file = Path(custom_ranges_dir, "custom_ranges_{}.txt".format(range))
    density_file = Path(custom_ranges_dir, "custom_ranges_{}.txt".format(range))

    symm_group_file = Path(custom_ranges_dir, "symm_groups_{}.txt".format(range))

    if mode == "cpu":
        mode_arg = "cpu_omp"
    elif mode == "gpu":
        mode_arg = "cuda"

    thresh_arg = ""
    if thresh:
        thresh_arg = "--skip -ct {}".format(thresh)

    cmd = "python3 {}/imp-sampcon/pyext/src/exhaust.py --sysname mtorc2 --path {} --rmfA {} --rmfB {} --scoreA {} --scoreB {} --selection {} --density {} --gnuplot --align --ambiguity {} -m {} {} -c 4 -g 2.0 -gp".format(
        Path.home(), analysis_dir, rmf_A_file, rmf_B_file, scores_A_file,
        scores_B_file, selection_file, density_file, symm_group_file, mode_arg, thresh_arg)

    print(cmd)
    os.system(cmd)

    print("Finished in {}s".format(time.time()-t0))


if __name__ == "__main__":
    analysis_dir = Path(sys.argv[1])
    sample_dir = Path(sys.argv[2])
    job_id = os.environ["JOB_ID"] if sys.argv[3] == "" else sys.argv[3]
    cluster = int(sys.argv[4])
    N_struct = int(sys.argv[5])
    filter = sys.argv[6]
    gsms_folder = sys.argv[7]
    range = sys.argv[8]
    copy = None if sys.argv[9] == "" else sys.argv[9]
    thresh = None if sys.argv[10] == "" else int(sys.argv[10])

    mode = "gpu"

    print("analysis_dir: {}".format(analysis_dir))
    print("sample_dir:   {}".format(sample_dir))
    print("job_id:       {}".format(job_id))
    print("cluster:      {}".format(cluster))
    print("N_struct:     {}".format(N_struct))
    print("filter:       {}".format(filter))
    print("gsms_folder:  {}".format(gsms_folder))
    print("range:        {}".format(range))
    print("copy:         {}".format(copy))
    print("thresh:       {}".format(thresh))
    print("mode:         {}".format(mode))

    sampcon_dir = Path(analysis_dir, "sampcon_{}".format(cluster))
    if not sampcon_dir.exists():
        sampcon_dir.mkdir()
        Path(sampcon_dir, "gsms").mkdir()

    job_dir = Path(sampcon_dir, job_id)
    job_dir.mkdir()

    gsms_dir = Path(sampcon_dir, "gsms/{}".format(gsms_folder))

    if not gsms_dir.exists():
        gsms_dir.mkdir()
        extract_models(
            output_dir=sample_dir,
            analysis_dir=analysis_dir,
            filter=filter,
            gsms_dir=gsms_dir,
            N=N_struct,
            cluster=cluster
        )
    else:
        print("gsms dir found")

    # # Need to copy over 3 files from an older sampcon run. The pre-computed distance matrix is used in clustering. Sampling precision statistics and Chi square grid statistics are for constructing plots.
    if copy:
        copy_files = list()
        copy_files.append("Distances_Matrix.data.npy")
        copy_files.append("mtorc2.Sampling_Precision_Stats.txt")
        copy_files.append("mtorc2.ChiSquare_Grid_Stats.txt")

        for file_name in copy_files:
            shutil.copy(Path(sampcon_dir, copy, file_name), Path(job_dir, file_name))

    os.chdir(job_dir)
    sampcon(
        analysis_dir=analysis_dir,
        gsms_dir=gsms_dir,
        mode=mode,
        thresh=thresh,
        range=range
    )