import sys, os
sys.path.append("/home/matthew/mtorc2/single_traj/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/single_traj/src")
from run_58 import run_58
from pathlib import Path


if __name__ == "__main__":
    print("Running modeling.py")
    tmp_dir = os.getenv('TMPDIR', '1')
    orig_dir = os.getenv('SGE_CWD_PATH', '1')
    job_id = os.getenv('JOB_ID', '1')
    print("tmp_dir: ", tmp_dir)
    print("orig_dir: ", orig_dir)
    print("job_id: ", job_id)

    xls = list()
    xls.append(("57.csv", 1, "DSS"))
    xls.append(("A.csv", 1, "DSS"))
    xls.append(("akt.csv", 1, "DSS"))
    max_dr = 25
    n_frames = 50000
    em_weight = 25
    print("max_dr: ", max_dr)
    print("n_frames: ", n_frames)
    print("em_weight: ", em_weight)

    output_dir = Path(tmp_dir, "output.{}".format(job_id), "output_0")
    print("output_dir: ", output_dir)

    density_dir = Path(Path.home(), "mtorc2/data/em/components")
    print("cp {}/* .".format(density_dir))
    os.system("cp {}/* .".format(density_dir))

    run_58(
        output_dir=output_dir,
        xls=xls,
        max_dr=max_dr,
        em_weight=em_weight,
        n_frames=n_frames,
    )

    print("mv output.{} {}".format(job_id, Path(orig_dir, "sample")))
    os.system("mv output.{} {}".format(job_id, Path(orig_dir, "sample")))
