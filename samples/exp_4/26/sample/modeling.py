import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from rigid_body_dimer import rigid_body_dimer


if __name__ == "__main__":
    xls = list()

    # Inter RB XLs: 10
    # MTOR intra RB: 2
    # RICTOR intra RB: 2
    xls.append(("inter_rb_dss.csv", 10, "DSS"))
    xls.append(("inter_rb_edc.csv", 10, "EDC"))
    xls.append(("mtor_intra_rb.csv", 2, "DSS"))
    xls.append(("RICTOR_intra_rb.csv", 2, "DSS"))

    t = rigid_body_dimer(
        xls=xls,
        n_steps=100000
    )
    print("Finished in {}s".format(t))