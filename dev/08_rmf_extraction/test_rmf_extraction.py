import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories
import shutil


if __name__ == "__main__":
    output_dir = Path("/wynton/home/sali/mhancock/mtorc2/samples/exp_14/125/sample/output.1894911")
    analysis_dir = Path("/wynton/group/sali/mhancock/mtorc2/samples/exp_14/125/analysis/1")
    cluster = 1
    gsms_dir = Path(Path.home(), "mtorc2/dev/08_rmf_extraction/gsms")
    N = 1000

    print("**Running extract_models***")
    t0 = time.time()

    out_dirs = [str(out_dir) for out_dir in output_dir.glob("*") if
                out_dir.is_dir()]

    # gsms_dir = Path(sampcon_dir, "gsms")
    hdbscan_A_file = Path(analysis_dir, "traj/A_cluster{}.csv".format(cluster))
    gsms_A_file = Path(gsms_dir, "A.csv".format(cluster, N))

    # Load module
    AT = AnalysisTrajectories(
        out_dirs,
        dir_name='output_',
        analysis_dir=str(Path(analysis_dir, "traj")),
        nproc=50
    )

    HA = AT.get_models_to_extract(str(gsms_A_file))

    rmf_A_file = Path(gsms_dir, "A.rmf3")
    scores_A_file = Path(gsms_dir, "A.txt")

    AT.extract_models_to_single_rmf(
        HA,
        str(rmf_A_file),
        str(output_dir),
        str(scores_A_file)
    )

    print("Finished in {}s".format(time.time()-t0))