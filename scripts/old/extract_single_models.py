import sys
import time
from pathlib import Path
sys.path.append("/home/matthew/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories


if __name__ == "__main__":
    i = 42
    run = 1
    exp_dir = "exp_5/{}".format(i)
    job_dir = Path(Path.home(), "mtorc2/single_traj_experiments", exp_dir)
    sample_dir = Path(job_dir, "sample")
    analysis_dir = Path(job_dir, "analysis/traj")
    out_dirs = [str(out_dir) for out_dir in sample_dir.glob("*") if out_dir.is_dir()]
    print(out_dirs)

    XLs_cutoffs = dict()
    XLs_cutoffs['DSS']=35.0
    # XLs_cutoffs['EDC']=35.0

    AT = AnalysisTrajectories(
        out_dirs,
        dir_name='output_',
        analysis_dir=str(analysis_dir),
        nproc=8,
        burn_in_fraction=.02,
        detect_equilibration=True,
        nskip=5
    )

    H = AT.get_models_to_extract(str(Path(analysis_dir, "scores_info_{}.csv".format(run))))
    AT.do_extract_models(
        gsms_info=H,
        filename="gsm_{}".format(run),
        gsms_dir=str(Path(job_dir, "analysis/gsms"))
    )