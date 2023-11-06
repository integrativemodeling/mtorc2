import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from rigid_body_dimer import rigid_body_dimer


if __name__ == "__main__":
    xls = list()
    xls.append(("inter_rb_dss.csv", 5, "DSS"))
    xls.append(("inter_rb_edc.csv", 5, "EDC"))

    t = rigid_body_dimer(
        xls=xls,
        n_steps=100000
    )
    print("Finished in {}s".format(t))