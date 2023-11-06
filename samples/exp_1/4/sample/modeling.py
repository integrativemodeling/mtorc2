import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from rigid_body_dimer import rigid_body_dimer


if __name__ == "__main__":
    t = rigid_body_dimer(
        w_rb_xls=2,
        n_steps=20000
    )
    print("Finished in {}s".format(t))