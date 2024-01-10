import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
import argparse
import multiprocessing

sys.path.append(str(Path(Path.home(), "PMI_analysis/pyext/src")))
from analysis_trajectories import AnalysisTrajectories
sys.path.append(str(Path(Path.home(), "mtorc2/src")))
import params


def analyze_trajectories(
        analysis_dir,
        out_dirs,
        cluster_on,
        equil,
        min_clust
):
    XLs_cutoffs = dict()
    XLs_cutoffs['DSS'] = 35.0
    XLs_cutoffs['EDC'] = 16.0

    print("EQUIL: {}".format(equil))
    if equil < 0:
        burn_in_fraction=.02
        detect_equilibration=True
    else:
        burn_in_fraction=equil
        detect_equilibration=False

    AT = AnalysisTrajectories(
        out_dirs,
        dir_name='output_',
        analysis_dir=str(Path(analysis_dir, "traj")),
        nproc=50,
        burn_in_fraction=burn_in_fraction,
        detect_equilibration=detect_equilibration,
        nskip=25
    )

    AT.set_analyze_XLs_restraint(
        XLs_cutoffs=XLs_cutoffs,
        ambiguous_XLs_restraint=True,
        get_nuisances=True
    )
    AT.set_analyze_Distance_restraint()
    AT.set_analyze_Excluded_volume_restraint()
    AT.set_analyze_EM_restraint()

    AT.read_stat_files()
    AT.write_models_info()
    AT.get_psi_stats()

    AT.hdbscan_clustering(
        cluster_on,
        min_cluster_size=min_clust
    )
    AT.summarize_XLs_info()


"""
AnalysisTrajectories density based clustering produces several clusters. The structures contained within each cluster are enumerated by an A and B csv file. In order to extract the structures for a given cluster, the A and B csv files must be corrected for the following reasons.

1) The column containing the rmf3 file for each structure is mangled. I am not sure why this is the case.

/scratch/9752918.110.member.q/output_109/rmfs//1.rmf3
/wynton/group/sali/mhancock/mtorc2/samples/exp_15/136/9752918/output_0/rmfs/

2) In order to combine multiple samples in the analysis pipeline, I add all output directories to a temporary head directory, but renumber all the output directories randomly. The output field for each entry therefore must be updated.
"""
def fix_rmf_file_pool(
        params_dict
):
    rmf_file_name = params_dict["rmf_file"].split("/")[-1]
    rmf_file_fixed = Path(params_dict["sample_dir"], params_dict["output_dir_name"], "rmfs", rmf_file_name)

    return rmf_file_fixed


def correct_csvs(
        analysis_dir,
        sample_dir
):
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    csv_files = list(Path(analysis_dir, "traj").glob("*_detailed.csv*"))

    print(analysis_dir)
    print(csv_files)
    for csv_file in csv_files:
        print(csv_file)
        cluster_df = pd.read_csv(csv_file)

        pool_params = list()
        rmf_files = cluster_df["rmf3_file"].tolist()
        output_dir_names = cluster_df["traj"].tolist()

        for i in range(len(rmf_files)):
            params_dict = dict()
            params_dict["rmf_file"] = rmf_files[i]
            params_dict["output_dir_name"] = output_dir_names[i]
            params_dict["sample_dir"] = sample_dir
            pool_params.append(params_dict)

        pool_results = pool_obj.imap(
            fix_rmf_file_pool,
            pool_params
        )

        cluster_df["rmf3_file"] = [str(rmf_file_fixed) for rmf_file_fixed in pool_results]

        comps_csv = str(csv_file.stem).split("_")

        half = comps_csv[2]
        cluster_name = comps_csv[3]

        cluster_df = cluster_df.rename(columns={'Unnamed: 0': ''})
        cluster_df.to_csv(Path(analysis_dir, "traj/{}_{}.csv".format(half, cluster_name)), index=False)


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis_dir", type=str, required=True)
    parser.add_argument("--sample_dir", type=str, required=True)
    parser.add_argument("--cluster_on", type=str, required=True)
    parser.add_argument("--equil", type=float, required=False)
    parser.add_argument("--min_clust", type=int, required=True)
    args = parser.parse_args()
    params.write_params(vars(args), Path(args.analysis_dir, "traj/params.txt"))

    analysis_dir = Path(args.analysis_dir)
    sample_dir = Path(args.sample_dir)
    cluster_on = args.cluster_on.split(",")

    out_dirs = [str(out_dir) for out_dir in sample_dir.glob("output_*")]
    # out_dirs = out_dirs[:10]

    t0 = time.time()
    analyze_trajectories(
        analysis_dir=analysis_dir,
        out_dirs=out_dirs,
        cluster_on=cluster_on,
        equil=args.equil,
        min_clust=args.min_clust
    )
    print("Analyze trajectories in {}s".format(time.time()-t0))

    t0 = time.time()
    correct_csvs(
        analysis_dir=analysis_dir,
        sample_dir=sample_dir
    )
    print("Corrected csvs in {}s".format(time.time()-t0))


