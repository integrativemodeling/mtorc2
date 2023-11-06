import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories


def analyze_trajectories(
        analysis_dir,
        out_dirs,
        cluster_on,
        equil,
        min_clust
):
    print("**Running analyze_trajectories***")
    print("analysis_dir: {}".format(analysis_dir))
    print("out_dirs:     {}".format(out_dirs))
    print("cluster_on:   {}".format(cluster_on))
    print("equil:        {}".format(equil))
    print("min_clust:    {}".format(min_clust))

    t0 = time.time()

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

    print("Finished in {}s".format(time.time()-t0))


"""
AnalysisTrajectories density based clustering produces several clusters. The structures contained within each cluster are enumerated by an A and B csv file. In order to extract the structures for a given cluster, the A and B csv files must be corrected for the following reasons.

1) The column containing the rmf3 file for each structure is mangled. I am not sure why this is the case. 

2) In order to combine multiple samples in the analysis pipeline, I add all output directories to a temporary head directory, but renumber all the output directories randomly. The output field for each entry therefore must be updated. 
"""
def correct_csvs(
        analysis_dir,
        output_dir
):
    print("**Running correct_csvs***")
    t0 = time.time()

    # csv_files = list(Path(analysis_dir, "traj").glob("selected*"))
    csv_files = list(Path(analysis_dir, "traj").glob("*_detailed.csv*"))

    for csv_file in csv_files:
        print(csv_file)
        cluster_df = pd.read_csv(csv_file)

        for i in range(len(cluster_df)):
            rmf_file = cluster_df.iloc[i, cluster_df.columns.get_loc("rmf3_file")]
            comps = rmf_file.split("/")

            new_traj_num = cluster_df.iloc[i, cluster_df.columns.get_loc("traj")]

            fixed_rmf_file = str(
                Path(output_dir, new_traj_num, comps[4], comps[6]))
            cluster_df.iloc[
                i, cluster_df.columns.get_loc("rmf3_file")] = fixed_rmf_file

        comps_csv = str(csv_file.stem).split("_")
        csv_sample_name = "{}_{}.csv".format(comps_csv[2],comps_csv[3])

        cluster_df = cluster_df.rename(columns={'Unnamed: 0': ''})
        cluster_df.to_csv(Path(analysis_dir, "traj", csv_sample_name), index=False)

    print("Finished in {}s".format(time.time()-t0))


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")

    analysis_dir = Path(sys.argv[1])
    sample_dir = Path(sys.argv[2])
    # os.system("mkdir {}".format(analysis_dir))
    # os.system("mkdir {}".format(Path(analysis_dir, "traj")))

    out_dirs = [str(out_dir) for out_dir in sample_dir.glob("output_*")]
    # out_dirs = out_dirs[:10]

    analyze_trajectories(
        analysis_dir=analysis_dir,
        out_dirs=out_dirs,
        cluster_on=sys.argv[3].split(","),
        equil=-1 if sys.argv[4] == "" else float(sys.argv[4]),
        min_clust=int(sys.argv[5])
    )

    correct_csvs(
        analysis_dir=analysis_dir,
        output_dir=sample_dir
    )


