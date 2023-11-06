import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from dimer_flex import dimer_flex


if __name__ == "__main__":

    t = dimer_flex(
        w_rb_xls=1,
        w_flex_xls=1,
        n_steps=20000
    )
    print("Finished in {}s".format(t))