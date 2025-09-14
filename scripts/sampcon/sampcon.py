import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
import shutil
import argparse

sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories
sys.path.append(str(Path(Path.home(), "mtorc2/src")))
import params


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

    out_dirs = [str(out_dir) for out_dir in output_dir.glob("*") if out_dir.is_dir()]

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

    print(HA.head())

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
        range_file,
        symm_file
):
    print("**Running sampcon***")
    t0 = time.time()

    rmf_A_file = Path(gsms_dir, "A.rmf3")
    rmf_B_file = Path(gsms_dir, "B.rmf3")
    scores_A_file = Path(gsms_dir, "A.txt")
    scores_B_file = Path(gsms_dir, "B.txt")

    # custom_ranges_dir = Path(Path.home(), "mtorc2/data/custom_ranges")
    # selection_file = Path(custom_ranges_dir, "custom_ranges_{}.txt".format(range))
    # density_file = Path(custom_ranges_dir, "custom_ranges_{}.txt".format(range))

    selection_file = range_file
    density_file = range_file

    # symm_group_file = Path(custom_ranges_dir, "symm_groups_{}.txt".format(range))

    if mode == "cpu":
        mode_arg = "cpu_omp"
    elif mode == "gpu":
        mode_arg = "cuda"

    thresh_arg = ""
    if thresh:
        thresh_arg = "--skip -ct {}".format(thresh)

    # cmd = "python3 {}/imp-sampcon/pyext/src/exhaust.py --sysname mtorc2 --path {} --rmfA {} --rmfB {} --scoreA {} --scoreB {} --selection {} --density {} --gnuplot --align --ambiguity {} -m {} {} -c 4 -g 2.0 -gp".format(
    #     Path.home(), analysis_dir, rmf_A_file, rmf_B_file, scores_A_file,
    #     scores_B_file, selection_file, density_file, symm_file, mode_arg, thresh_arg)

    cmd = "python3 {}/imp-sampcon/pyext/src/exhaust.py --sysname mtorc2 --path {} --rmfA {} --rmfB {} --scoreA {} --scoreB {} --selection {} --density {} --gnuplot --align -m {} {} -c 4 -g 2.0 -gp".format(
        Path.home(), analysis_dir, rmf_A_file, rmf_B_file, scores_A_file,
        scores_B_file, selection_file, density_file, mode_arg, thresh_arg)

    if symm_file:
        cmd += " --ambiguity {}".format(symm_file)

    print(cmd)
    os.system(cmd)

    print("Finished in {}s".format(time.time()-t0))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis_dir", type=str, required=True)
    parser.add_argument("--sample_dir", type=str, required=True)
    parser.add_argument("--job_id", type=str, required=True)
    parser.add_argument("--score_cluster", type=int, required=True)
    parser.add_argument("--n_struct", type=int, required=True)
    parser.add_argument("--filter", type=str, required=True)
    parser.add_argument("--gsms_dir_name", type=str, required=True)
    parser.add_argument("--range_file", type=str, required=True)
    parser.add_argument("--symm_file", type=str, required=False)
    parser.add_argument("--copy", type=int, required=False)
    parser.add_argument("--thresh", type=int, required=False)
    parser.add_argument("--extract_only", action="store_true")
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    sample_dir = Path(args.sample_dir)

    range_file = Path(args.range_file)

    print(args.copy)
    print(args.thresh)

    mode = "cpu"

    sampcon_dir = Path(analysis_dir, "sampcon_{}".format(args.score_cluster))
    if not sampcon_dir.exists():
        sampcon_dir.mkdir()
        Path(sampcon_dir, "gsms").mkdir()

    job_dir = Path(sampcon_dir, args.job_id)
    gsms_dir = Path(sampcon_dir, "gsms/{}".format(args.gsms_dir_name))

    if not gsms_dir.exists():
        print("Creating good scoring models directory")
        gsms_dir.mkdir()
        extract_models(
            output_dir=sample_dir,
            analysis_dir=analysis_dir,
            filter=args.filter,
            gsms_dir=gsms_dir,
            N=args.n_struct,
            cluster=args.score_cluster
        )
    else:
        print("Good scoring models directory found")

    if not args.extract_only:
        # # #
        if job_dir.exists():
            shutil.rmtree(job_dir)
        # if gsms_dir.exists():
        #     shutil.rmtree(gsms_dir)
        # # #

        # # Raise an error if the job_dir already exists.
        job_dir.mkdir(exist_ok=False)
        params.write_params(vars(args), Path(job_dir, "params.txt"))

        # Need to copy over 3 files from an older sampcon run. The pre-computed distance matrix is used in clustering. Sampling precision statistics and Chi square grid statistics are for constructing plots.
        if type(args.copy) == int:
            copy_from_dir = Path(sampcon_dir, str(args.copy))
            orig_files = [Path(copy_from_dir, "Distances_Matrix.data.npy"), Path(copy_from_dir, "mtorc2.Sampling_Precision_Stats.txt"), Path(copy_from_dir, "mtorc2.ChiSquare_Grid_Stats.txt")]

            for orig_file in orig_files:
                new_file = Path(job_dir, orig_file.name)
                print(orig_file, new_file)
                shutil.copy(orig_file, new_file)

        os.chdir(job_dir)
        sampcon(
            analysis_dir=analysis_dir,
            gsms_dir=gsms_dir,
            mode=mode,
            thresh=args.thresh,
            range_file=range_file,
            symm_file=args.symm_file
        )
    else:
        print("Extracting GSMs only")