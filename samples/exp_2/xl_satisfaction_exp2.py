import pandas as pd
from pathlib import Path
import multiprocessing


def get_n_satisfied(XL_data, inter_vs_intra=False, ambig=False):
    n_satisfied = 0
    cutoff = 35
    if ambig:
        n_unique_entries = int(len(XL_data) / 4)
        for i in range(n_unique_entries):
            start=i*4
            end=(i+1)*4
            if XL_data[start:end].min() < cutoff:
                n_satisfied = n_satisfied + 1
    else:
        for i in range(len(XL_data)):
            if XL_data[i] < cutoff:
                n_satisfied = n_satisfied + 1

    return n_satisfied


def run_analysis(id):
    exp_num = "exp_2/{}".format(id)
    ambig = True
    inter_vs_intra = False
    file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/traj/XLs_dist_info_0.csv")
    df = pd.read_csv(file)
    df = df.drop(columns=["Unnamed: 0"])
    # print(df.columns)
    df.head()
    print(len(df))

    n_sats = list()
    for i in range(len(df)):
        # if i % 10000 == 0:
        #     print(i)
        XL_data = df.iloc[i,1:]
        n_sat = get_n_satisfied(XL_data, inter_vs_intra=inter_vs_intra, ambig=ambig)
        n_sats.append(n_sat)

    df["n_sat"] = n_sats
    # df.head()

    max = df.n_sat.max()
    n_frames = len(df[df.n_sat == df.n_sat.max()])
    print("MAX: ", max)
    print("N FRAMES: ", n_frames)
    print(df[df.n_sat == df.n_sat.max()])

    return id, max, n_frames


if __name__ == "__main__":
    # Compute all the scores in parallel.
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    score_results = pool_obj.imap(
        run_analysis,
        list(range(10,15))
    )

    r_factors_dict = dict()
    for id, max, n_frames in score_results:
        print(id, max, n_frames)
    pool_obj.close()

    # Write all decoy scores to a pickle file.
    # with open(scores_file, 'wb') as file:
    #     pickle.dump(r_factors_dict, file)





