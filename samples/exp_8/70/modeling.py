import sys, os
sys.path.append("/home/matthew/mtorc2/single_traj/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/single_traj/src")
from run_65 import run_65
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
    xls.append(("B.csv", 1, "DSS"))
    xls.append(("akt.csv", 1, "DSS"))
    n_frames = 100000
    print("n_frames: ", n_frames)

    output_dir = Path(tmp_dir, "output.{}".format(job_id), "output_0")
    print("output_dir: ", output_dir)

    density_dir = Path(Path.home(), "mtorc2/data/em/components_1")
    print("cp {}/* .".format(density_dir))
    os.system("cp {}/* .".format(density_dir))

    run_65(
        output_dir=output_dir,
        xls=xls,
        n_frames=n_frames,
        em=False
    )

    print("mv output.{} {}".format(job_id, Path(orig_dir, "sample")))
    os.system("mv output.{} {}".format(job_id, Path(orig_dir, "sample")))
