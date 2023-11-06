import sys
sys.path.append("/home/matthew/mtorc2/single_traj_experiments/src")
from dimer_flex import dimer_flex


if __name__ == "__main__":
    xls = list()
    xls.append(("inter_rb_dss.csv", 1, "DSS"))
    xls.append(("inter_rb_edc.csv", 1, "EDC"))
    xls.append(("inter_flex_dss.csv", 1, "DSS"))
    xls.append(("inter_flex_edc.csv", 1, "EDC"))
    xls.append(("mtor_intra_rb.csv", 2, "DSS"))
    xls.append(("mtor_intra_flex.csv", 1, "DSS"))

    t = dimer_flex(
        xls=xls,
        n_steps=100000
    )
    print("Finished in {}s".format(t))