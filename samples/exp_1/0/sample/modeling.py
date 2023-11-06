import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from rigid_body_monomer import rigid_body_monomer
from pathlib import Path


if __name__ == "__main__":
    t = rigid_body_monomer(
        w_rb_xls=1,
        n_steps=20000
    )
    print("Finished in {}s".format(t))