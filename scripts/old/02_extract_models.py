import sys, os
from pathlib import Path
sys.path.append("/home/matthew/PMI_analysis/pyext/src")
sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories


if __name__ == "__main__":
    c = 3 # Cluster we wish to extract
    N = 1000

    exp_info = (8, 65, 90001, ["57", "A", "akt"])

    output_dir = Path(Path.home(), "mtorc2/exp_{}/{}/sample/output.{}".format(exp_info[0], exp_info[1], exp_info[2]))
    analysis_dir = Path(Path.home(), "mtorc2/exp_{}/{}/analysis/{}".format(exp_info[0], exp_info[1], exp_info[2]))

    out_dirs = [str(out_dir) for out_dir in output_dir.glob("*") if out_dir.is_dir()]
    print(out_dirs)
    cluster_dir = Path(analysis_dir, "cluster_{}".format(c))
    os.system("mkdir {}".format(cluster_dir))

    # Load module
    AT = AnalysisTrajectories(
        out_dirs,
        dir_name='output_',
        analysis_dir=str(analysis_dir),
        nproc=8,
        burn_in_fraction=.02,
        detect_equilibration=True,
        nskip=100
    )

    # Point to the selected_models file
    A_file = Path(analysis_dir, "A_cluster{}_{}.csv".format(c, N))
    B_file = Path(analysis_dir, "B_cluster{}_{}.csv".format(c, N))
    HA = AT.get_models_to_extract(str(A_file))
    HB = AT.get_models_to_extract(str(B_file))

    rmf_A_file = Path(cluster_dir, "A_{}.rmf3".format(N))
    rmf_B_file = Path(cluster_dir, "B_{}.rmf3".format(N))

    scores_A_file = Path(cluster_dir, "A_score_{}.txt".format(N))
    scores_B_file = Path(cluster_dir, "B_score_{}.txt".format(N))

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

