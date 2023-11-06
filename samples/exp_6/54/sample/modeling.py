from pathlib import Path
import sys, os
sys.path.append("/home/matthew/mtorc2/single_traj/src")
sys.path.append("/wynton/home/sali/mhancock/mtorc2/single_traj/src")
from run_54 import run_54


if __name__ == "__main__":
    xls = list()
    xls.append(("45.csv", 1, "DSS"))
    max_dr = None
    n_frames = 25000

    if len(sys.argv) > 1:
        home_dir = Path(str(sys.argv[1]))
    else:
        home_dir = Path.home()

    exp_dir = Path(home_dir, "mtorc2/single_traj", Path.cwd().parents[1].name)
    output_dir = Path(exp_dir, Path.cwd().parents[0].name, "sample/output_0")
    print("OUT DIR:", output_dir)

    run_54(
        output_dir=output_dir,
        xls=xls,
        max_dr=max_dr,
        n_frames=n_frames,
    )

    if len(sys.argv) > 1:
        print("cp -r {} .".format(output_dir))
        os.system("cp -r {} .".format(output_dir))
