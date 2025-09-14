from pathlib import Path
import pandas as pd
from pathlib import Path

import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import RMF
import matplotlib.pyplot as plt
import multiprocessing


def get_n_clash_pool(
    pool_input
):
    rmf_file, frame_id = pool_input
    print(rmf_file, frame_id)

    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(frame_id))

    pids = IMP.atom.Selection(h, resolution=10).get_selected_particle_indexes()

    n_clash = 0
    n_total = len(pids) * (len(pids)-1) / 2
    for pid_0 in pids:
        xyzr_0 = IMP.core.XYZR(m, pid_0)
        for pid_1 in pids:
            xyzr_1 = IMP.core.XYZR(m, pid_1)

            if pid_0 != pid_1 and IMP.core.get_distance(xyzr_0, xyzr_1) < xyzr_0.get_radius() + xyzr_1.get_radius():
                n_clash += 1


    return rmf_file, frame_id, n_clash, n_total


if __name__ == "__main__" :
    gsms_dir = Path("/wynton/group/sali/mhancock/mtorc2/analysis/136/2/sampcon_-1/gsms/1000R")

    pool_inputs = list()
    for rmf_file in [Path(gsms_dir, "A.rmf3"), Path(gsms_dir, "B.rmf3")]:
        for i in range(1000):
            pool_inputs.append((rmf_file, i))

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())

    score_results = pool_obj.imap(
        get_n_clash_pool,
        pool_inputs
    )

    clash_df = pd.DataFrame(columns=["rmf_file", "frame", "n_clash", "n_total"])
    for rmf_file, frame_id, n_clash, n_total in score_results:
        clash_df.loc[len(clash_df)] = [rmf_file, frame_id, n_clash, n_total]

    pool_obj.close()

    clash_df.to_csv(Path(Path.home(), "mtorc2/data/models/submission_3/136_2_1/ev_clash.csv"))