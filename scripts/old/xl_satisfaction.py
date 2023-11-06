import pandas as pd
from pathlib import Path
import multiprocessing
import IMP
import IMP.atom
import IMP.rmf
import RMF


def tuple_to_string(entry):
    return entry[0] + "_" + str(entry[1]) + "_" + entry[2] + "_" + str(entry[3])


def string_to_tuple(entry):
    comps = entry.split("_")
    return comps[0], int(comps[1]), comps[2], int(comps[3])


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


# Return true if the model contains the following residue at any resolution.
def contains_residue(m, res, prot):
    particle_names = list()
    for pid in m.get_particle_indexes():
        if m.get_particle_name(pid) in ["MTOR", "RICTOR", "MLST8", "MSIN1", "CRIM", "RBD", "AKT1"]:
            mol = m.get_particle_name(pid)

        if IMP.atom.Fragment.get_is_setup(m, pid):
            frag = IMP.atom.Fragment(m, pid)
            if frag.get_contains_residue(res) and mol == prot:
                particle_names.append(m.get_particle_name(pid))

    if len(particle_names) > 0:
        return True
    else:
        return False


# XL sat produces 2 pandas dataframe. The metadata dataframe (all_xls_df) contains information about the best scoring frame. Important is whether the best scoring model contains a given XL ('contains'), whether that XL was satisfied ('satisfied'), and the copy numbers for the minimum XL distance ('copy1', 'copy2'). The XL distance dataframe contains the minimum XL distance for each XL that the model was optimized against.
def generate_xl_sat_dfs(exp_num):
    # First get the best scoring frame number.
    scores_file = Path(Path.home(), "mtorc2/single_traj_experiments", exp_num, "analysis/traj/scores_info_0.csv")
    scores_df = pd.read_csv(scores_file)
    best_scoring_frame = scores_df[scores_df["Total_Score"] == scores_df["Total_Score"].min()].index[0]

    # Read in the raw XL distance info from PMI analysis.
    raw_xl_dist_file = Path(Path.home(), "mtorc2/single_traj_experiments", exp_num, "analysis/traj/XLs_dist_info_0.csv")
    raw_xl_df = pd.read_csv(raw_xl_dist_file)
    raw_xl_df = raw_xl_df.drop(columns=["Unnamed: 0"])

    # Read in the XL info.
    all_xls_file = Path(Path.home(), "mtorc2/data/xlms/monomer.xl_satisfaction.csv")
    all_xls_df = pd.read_csv(all_xls_file)
    all_xls_df = all_xls_df.drop(columns=["Unnamed: 0"])
    all_xls_df.satisfied = False
    all_xls_df.contains = False
    all_xls_df.head()

    # Add fields for the copy id for the best scoring structure.
    all_xls_df["copy1"] = 0
    all_xls_df["copy2"] = 0

    # Iterate through all XLs and find XLs that were included in modeling. If a XL is included in the model; determine what are the column IDs of the ambiguous XL distances as recorded in the raw XL distance df (matched_cols_dict). An example entry is [1,2,3,4] Additionally, record the copy numbers for all matched columns (copy_orders_dict).  An example entry is [(0,0), (1,0), (1,1), (0,1)].
    col_strs = list()
    copy_orders_dict = dict()
    matched_cols_dict = dict()
    for i in range(len(all_xls_df)):
        prot1, res1, prot2, res2 = tuple(all_xls_df.iloc[i, 0:4])
        entry = (prot1, res1, prot2, res2)

        copy_orders_dict[entry] = list()
        matched_cols = get_columns(raw_xl_df, prot1, res1, prot2, res2)
        if len(matched_cols) > 0:
            col_strs.append(tuple_to_string(entry))
            matched_cols_dict[entry] = matched_cols
            for col in matched_cols:
                copy_1 = raw_xl_df.columns[col].split('|')[3].split(".")[1]
                copy_2 = raw_xl_df.columns[col].split('|')[5].split(".")[1]
                copy_orders_dict[entry].append((copy_1, copy_2))

    # Use both dictionaries to create a new dataframe that removes the ambiguity from the raw XL distance dataframe (xl_sat_df). All entries in the new dataframe are unique and represent the minimum XL distance for a given entry.
    xl_sat_dict = dict()
    for entry in matched_cols_dict.keys():
        xl_sat_dict[tuple_to_string(entry)] = list()

    best_scoring_dict = dict()
    for i in range(len(raw_xl_df)):
        for col_str in col_strs:
            entry = string_to_tuple(col_str)
            col_ids = matched_cols_dict[entry]
            copy_orders = copy_orders_dict[entry]
            min_xl_dist = raw_xl_df.iloc[i, col_ids[0]]
            min_copy_id = ('0','0')
            for j in range(4):
                xl_dist = raw_xl_df.iloc[i, col_ids[j]]
                if xl_dist < min_xl_dist:
                    min_xl_dist = xl_dist
                    min_copy_id = copy_orders[j]

            xl_sat_dict[col_str].append(min_xl_dist)

            if i == best_scoring_frame:
                best_scoring_dict[entry] = (min_xl_dist, min_copy_id)

    xl_sat_df = pd.DataFrame(xl_sat_dict)
    xl_sat_file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/XLs_sat.csv")
    xl_sat_df.to_csv(xl_sat_file)

    # Update the dataframe containing the metadata for all XLs and save.
    rmf_file = Path(Path.home(), "Documents/mtorc2/single_traj_experiments/exp_5/37/sample/output_0/rmfs/0.rmf3")
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, 0)
    for entry in best_scoring_dict.keys():
        for i in range(len(all_xls_df)):
            # if entry == tuple(all_xls_df.iloc[i, 0:4]):
            prot1, res1, prot2, res2 = all_xls_df.iloc[i, 0:4]
            contains1 = contains_residue(m, prot1, res1)
            contains2 = contains_residue(m, prot2, res2)
            if contains1 and contains2:
                all_xls_df.iloc[i, all_xls_df.columns.get_loc("contains")] = True

                xl_dist, copy_id = best_scoring_dict[entry]
                if xl_dist < 35:
                    all_xls_df.iloc[i, all_xls_df.columns.get_loc("satisfied")] = True

                all_xls_df.iloc[i, all_xls_df.columns.get_loc("copy1")] = int(copy_id[0])
                all_xls_df.iloc[i, all_xls_df.columns.get_loc("copy2")] = int(copy_id[1])

    best_frame_file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num),"analysis/best_frame.csv")
    all_xls_df.to_csv(best_frame_file)

    return None

            # def get_n_satisfied(exp_num):
#     file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/traj/XLs_dist_info_0.csv")
#     df = pd.read_csv(file)
#     df = df.drop(columns=["Unnamed: 0"])
#     df.head()
#
#     init_file = Path(Path.home(), "mtorc2/data/xlms/monomer.xl_satisfaction.csv")
#     init_df = pd.read_csv(init_file)
#     init_df = init_df.drop(columns=["Unnamed: 0"])
#     init_df.satisfied = False
#     init_df.head()
#
#     cols = list()
#     for i in range(len(init_df)):
#         prot1 = init_df.iloc[i, init_df.columns.get_loc("prot1")]
#         prot2 = init_df.iloc[i, init_df.columns.get_loc("prot2")]
#         res1 = int(init_df.iloc[i, init_df.columns.get_loc("res1")])
#         res2 = int(init_df.iloc[i, init_df.columns.get_loc("res2")])
#
#         matched_cols = get_columns(df, prot1, res1, prot2, res2)
#         if len(matched_cols) > 0:
#             cols.append(tuple_to_string((prot1, res1, prot2, res2)))
#
#     cols_dict = dict()
#     col_ids_dict = dict()
#     for col in cols:
#         prot1, res1, prot2, res2 = string_to_tuple(col)
#         col_ids_dict[col] = get_columns(df, prot1, res1, prot2, res2)
#         cols_dict[col] = list()
#
#     for i in range(len(df)):
#         for col in cols:
#             col_ids = col_ids_dict[col]
#             min_xl_dist = df.iloc[i, col_ids[0]]
#             for col_id in col_ids:
#                 xl_dist = df.iloc[i, col_id]
#                 if xl_dist < min_xl_dist:
#                     min_xl_dist = xl_dist
#
#             cols_dict[col].append(min_xl_dist)
#
#     xl_sat_df = pd.DataFrame(cols_dict)
#
#     # intra_sats, inter_sats = list(), list()
#     # for i in range(len(xl_sat_df)):
#     #     n_intra_sat, n_inter_sat = 0, 0
#     #     for j in range(len(xl_sat_df.columns)):
#     #         prot1, res1, prot2, res2 = string_to_tuple(xl_sat_df.columns[j])
#     #         if xl_sat_df.iloc[i, j] < 35:
#     #             if prot1 == prot2:
#     #                 n_intra_sat = n_intra_sat+1
#     #             else:
#     #                 n_inter_sat = n_inter_sat+1
#     #
#     #     intra_sats.append(n_intra_sat)
#     #     inter_sats.append(n_inter_sat)
#     #
#     # xl_sat_df["inter_xl_sat"] = inter_sats
#     # xl_sat_df["intra_xl_sat"] = intra_sats
#
#     xl_sat_df.head()
#
#     xl_sat_file = Path(Path.home(), "mtorc2/single_traj_experiments", str(exp_num) ,"analysis/XLs_sat.csv")
#     xl_sat_df.to_csv(xl_sat_file)

    # inter_max = xl_sat_df.inter_xl_sat.max()
    # intra_max = xl_sat_df.intra_xl_sat.max()
    # n_max_inter = len(xl_sat_df[xl_sat_df.inter_xl_sat == inter_max])
    # n_max_intra = len(xl_sat_df[xl_sat_df.intra_xl_sat == intra_max])
    #
    # return exp_num, inter_max, n_max_inter, intra_max, n_max_intra

    # return None


if __name__ == "__main__":
    # Compute all the scores in parallel.
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    score_results = pool_obj.imap(
        generate_xl_sat_dfs,
        ["exp_5/" + str(i) for i in range(37, 38)]
    )

    r_factors_dict = dict()
    # for id, inter_max, n_max_inter, intra_max, n_max_intra in score_results:
    #     print(id, inter_max, n_max_inter, intra_max, n_max_intra)
    for result in score_results:
        continue

    pool_obj.close()

    # Write all decoy scores to a pickle file.
    # with open(scores_file, 'wb') as file:
    #     pickle.dump(r_factors_dict, file)





