from pathlib import Path
import pandas as pd
import random


def get_random_entries(
        df,
        N
):
    ids = list()
    for i in range(N):
        ids.append(int(random.random() * len(df)))

    sample_df = df.iloc[ids]
    sample_df.reset_index(
        inplace=True,
        drop=True
    )

    return sample_df


if __name__ == "__main__":
    exp_info = (8, 65, 90001, ["57", "A", "akt"])
    analysis_dir = Path(Path.home(), "mtorc2/exp_{}/{}/analysis/{}".format(exp_info[0], exp_info[1], exp_info[2]))

    # csv_files = analysis_dir.glob("selected_models_*random*")
    csv_files = list()
    csv_files.append("selected_models_A_cluster3_detailed.csv")
    csv_files.append("selected_models_B_cluster3_detailed.csv")
    csv_files = [Path(analysis_dir, csv_file) for csv_file in csv_files]
    n_structures = 1000

    sample_dir = Path(Path.home(), "mtorc2/exp_{}/{}/sample/output.{}".format(exp_info[0], exp_info[1], exp_info[2]))

    for csv_file in csv_files:
        # print(csv_file)
        cluster_df = pd.read_csv(csv_file)

        sample_df = get_random_entries(
            cluster_df,
            n_structures
        )

        for i in range(len(sample_df)):
            rmf_file = sample_df.iloc[i, sample_df.columns.get_loc("rmf3_file")]
            comps = rmf_file.split("/")
            # print(comps)
            fixed_rmf_file = str(Path(sample_dir, comps[3], comps[4], comps[6]))
            # print(fixed_rmf_file)
            sample_df.iloc[i, sample_df.columns.get_loc("rmf3_file")] = fixed_rmf_file

        comps_csv = str(csv_file.stem).split("_")

        csv_sample_name = comps_csv[2]+"_"+comps_csv[3]+"_{}".format(n_structures)+".csv"
        print(csv_sample_name)

        fixed_csv_file = Path(analysis_dir, csv_sample_name)

        sample_df = sample_df.rename(columns={'Unnamed: 0': ''})

        print(sample_df.head())

        sample_df.to_csv(fixed_csv_file, index=False)