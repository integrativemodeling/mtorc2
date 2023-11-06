from xl_satisfaction_cluster import check_ambiguity
from xl_satisfaction_cluster import check_satisfaction
from xl_satisfaction_cluster import draw_xls
import IMP
import IMP.atom
import IMP.rmf
import RMF

from pathlib import Path
import pandas as pd

if __name__ == "__main__":
    # rmf_file = Path(cluster_dir, "cluster_center_model.rmf3")
    rmf_file = Path("/wynton/home/sali/mhancock/mtorc2/exp_13/118/sample/test/output_0/initial.0.rmf3")
    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, RMF.FrameID(0))

    xl_file = Path(Path.home(), "mtorc2/data/xlms/all.csv")
    xl_df = pd.read_csv(xl_file)
    xl_df = xl_df.drop(columns=["Unnamed: 0"])

    xl_df["sat"] = 0
    xl_df["ambig"] = 0
    xl_df["frame"] = -1

    xl_df = check_satisfaction(
        xl_df=xl_df,
        m=m, h=h
    )

    xl_df = check_ambiguity(
        xl_df=xl_df,
        m=m, h=h
    )

    pb_file = Path(Path.home(), "mtorc2/data/xlms/visualization/80.pb")
    draw_xls(
        xl_df=xl_df,
        m=m, h=h,
        pb_file=pb_file
    )

