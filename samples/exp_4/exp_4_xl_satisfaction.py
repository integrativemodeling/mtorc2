import pandas as pd
from pathlib import Path
import multiprocessing


def get_columns(df, prot_1, res_1, prot_2, res_2):
    matched_col_ids = list()
    for i in range(1, len(df.columns)):
        curr_col = df.columns[i]
        curr_entry = curr_col.split('|')
        curr_prot_1 = curr_entry[3].split('.')[0]
        curr_prot_2 = curr_entry[5].split('.')[0]
        curr_res_1 = int(curr_entry[4])
        curr_res_2 = int(curr_entry[6])

        if prot_1 == curr_prot_1 and prot_2 == curr_prot_2 and res_1 == curr_res_1 and res_2 == curr_res_2:
            matched_col_ids.append(i)
        elif prot_2 == curr_prot_1 and prot_1 == curr_prot_2 and res_2 == curr_res_1 and res_1 == curr_res_2:
            matched_col_ids.append(i)

    return matched_col_ids


def get_n_satisfied(id):
    exp_num = "exp_3/{}".format(id)
    file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/traj/XLs_dist_info_0.csv")
    df = pd.read_csv(file)
    df = df.drop(columns=["Unnamed: 0"])
    print(len(df.columns))
    df.head()

    init_file = Path(Path.home(), "mtorc2/utility/initial_model/data/init_xl_satisfaction.csv")
    init_df = pd.read_csv(init_file)
    init_df = init_df.drop(columns=["Unnamed: 0"])
    init_df.satisfied = False
    init_df.head()

    cols = list()
    for i in range(len(init_df)):
        prot1 = init_df.iloc[i, init_df.columns.get_loc("prot1")]
        prot2 = init_df.iloc[i, init_df.columns.get_loc("prot2")]
        res1 = int(init_df.iloc[i, init_df.columns.get_loc("res1")])
        res2 = int(init_df.iloc[i, init_df.columns.get_loc("res2")])

        matched_cols = get_columns(df, prot1, res1, prot2, res2)
        if len(matched_cols) > 0:
            cols.append((prot1, res1, prot2, res2))

    cols_dict = dict()
    for col in cols:
        cols_dict[col] = list()

    col_ids_dict = dict()
    for col in cols:
        prot1, res1, prot2, res2 = col
        col_ids_dict[col] = get_columns(df, prot1, res1, prot2, res2)

    for i in range(len(df)):
        # if i % 10000 == 0:
        #     print(i)
        for col in cols:
            prot1, res1, prot2, res2 = col
            col_ids = col_ids_dict[col]
            #         print(col_ids)
            min_xl_dist = df.iloc[i, col_ids[0]]
            for col_id in col_ids:
                xl_dist = df.iloc[i, col_id]
                #             print(xl_dist)
                if xl_dist < min_xl_dist:
                    min_xl_dist = xl_dist

            cols_dict[col].append(min_xl_dist)

    xl_sat_df = pd.DataFrame(cols_dict)
    xl_sat_df.head()

    intra_sats, inter_sats = list(), list()
    for i in range(len(xl_sat_df)):
        # if i % 10000 == 0:
        #     print(i)

        n_intra_sat, n_inter_sat = 0, 0
        for j in range(len(xl_sat_df.columns)):
            prot1, res1, prot2, res2 = xl_sat_df.columns[j]
            if xl_sat_df.iloc[i, j] < 35:
                if prot1 == prot2:
                    n_intra_sat = n_intra_sat+1
                else:
                    n_inter_sat = n_inter_sat+1

        intra_sats.append(n_intra_sat)
        inter_sats.append(n_inter_sat)

    xl_sat_df["inter_xl_sat"] = inter_sats
    xl_sat_df["intra_xl_sat"] = intra_sats

    xl_sat_file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/XLs_sat.csv")
    xl_sat_df.to_csv(xl_sat_file)

    inter_max = xl_sat_df.inter_xl_sat.max()
    intra_max = xl_sat_df.intra_xl_sat.max()
    n_max_inter = len(xl_sat_df[xl_sat_df.inter_xl_sat == inter_max])
    n_max_intra = len(xl_sat_df[xl_sat_df.intra_xl_sat == intra_max])

    return id, inter_max, n_max_inter, intra_max, n_max_intra


if __name__ == "__main__":
    # Compute all the scores in parallel.
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    score_results = pool_obj.imap(
        get_n_satisfied,
        list(range(15,25))
    )

    r_factors_dict = dict()
    for id, inter_max, n_max_inter, intra_max, n_max_intra in score_results:
        print(id, inter_max, n_max_inter, intra_max, n_max_intra)
    pool_obj.close()

    # Write all decoy scores to a pickle file.
    # with open(scores_file, 'wb') as file:
    #     pickle.dump(r_factors_dict, file)





