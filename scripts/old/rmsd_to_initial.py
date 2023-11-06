import pandas as pd
from pathlib import Path
import multiprocessing
import IMP
import IMP.rmf
import IMP.atom
import RMF


def get_em_rmsds(exp_id):
    score_file = Path(Path.home(), "mtorc2/single_traj_experiments/exp_3", str(exp_id), "analysis/traj/scores_info_0.csv")
    save_file = Path(Path.home(), "mtorc2/single_traj_experiments/exp_3", str(exp_id), "analysis/em_rmsd.csv")
    score_df = pd.read_csv(score_file)
    score_df["XL_sum"] = score_df["XLs_DSS"] + score_df["XLs_EDC"]

    best_score_df = score_df.nsmallest(100, columns=["XL_sum"])
    em_rmsds = list()
    for index in list(best_score_df.index):
        print(index)
        em_rmsds.append(get_em_rmsd(exp_id, index))

    best_score_df["em_rmsd"] = em_rmsds
    best_score_df.to_csv(save_file)

    mean_rmsd = best_score_df["em_rmsd"].mean()
    return exp_id, mean_rmsd


def get_em_rmsd(exp_id, frame_num):
    file_1 = Path(Path.home(), "mtorc2/single_traj_experiments/exp_3", str(exp_id), "sample/output_0/initial.0.rmf3")
    file_2 = Path(Path.home(), "mtorc2/single_traj_experiments/exp_3", str(exp_id), "sample/output_0/rmfs/0.rmf3")

    ref_fh = RMF.open_rmf_file_read_only(str(file_1))
    sample_fh = RMF.open_rmf_file_read_only(str(file_2))
    ref_m = IMP.Model()
    sample_m = IMP.Model()
    ref_h = IMP.rmf.create_hierarchies(ref_fh, ref_m)[0]
    sample_h = IMP.rmf.create_hierarchies(sample_fh, sample_m)[0]
    IMP.rmf.load_frame(sample_fh, frame_num)

    ref_ps = IMP.atom.Selection(ref_h, molecules=['MTOR', 'MLST8', 'RICTOR'], resolution=1, copy_index=0).get_selected_particles()
    sample_ps = IMP.atom.Selection(sample_h, molecules=['MTOR', 'MLST8', 'RICTOR'], resolution=1, copy_index=0).get_selected_particles()

    t = IMP.atom.get_transformation_aligning_first_to_second(sample_ps, ref_ps)
    # IMP.atom.transform(sample_h, t) # transform the hierarchy in place
    rmsd = IMP.atom.get_rmsd_transforming_first(t, sample_ps, ref_ps)
    return rmsd


if __name__ == "__main__":
    # Compute all the scores in parallel.
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    score_results = pool_obj.imap(
        get_em_rmsds,
        list(range(15,25))
    )

    r_factors_dict = dict()
    for exp_id, mean_rmsd in score_results:
        print(exp_id, mean_rmsd)
    # pool_obj.close()

    # Write all decoy scores to a pickle file.
    # with open(scores_file, 'wb') as file:
    #     pickle.dump(r_factors_dict, file)





