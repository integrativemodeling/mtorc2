import sys
import os
import time
from pathlib import Path
import pandas as pd
import random
sys.path.append("/wynton/home/sali/mhancock/PMI_analysis/pyext/src")
from analysis_trajectories import AnalysisTrajectories
import shutil
import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.rmf
import RMF


def get_rmsd(
        rmf_1,
        frame_1,
        rmf_2,
        frame_2
):

    fh_1 = RMF.open_rmf_file_read_only(str(rmf_1))
    m_1 = IMP.Model()
    h_1 = IMP.rmf.create_hierarchies(fh_1, m_1)[0]
    IMP.rmf.load_frame(fh_1, RMF.FrameID(frame_1))

    fh_2 = RMF.open_rmf_file_read_only(str(rmf_2))
    m_2 = IMP.Model()
    h_2 = IMP.rmf.create_hierarchies(fh_2, m_2)[0]
    IMP.rmf.load_frame(fh_2, RMF.FrameID(frame_2))

    xyzs_1 = IMP.core.XYZs(IMP.core.get_leaves(h_1))
    xyzs_2 = IMP.core.XYZs(IMP.core.get_leaves(h_2))
    rmsd = IMP.atom.get_rmsd(xyzs_1, xyzs_2)

    print(rmsd)


if __name__ == "__main__":
    # output_dir = Path(Path.home(), "mtorc2/exp_10/94/sample/output.114524")
    # out_dirs = [Path(output_dir, "output_0")]
    # analysis_dir = Path(Path.home(), "mtorc2/exp_10/94/analysis/176064")
    # AT = AnalysisTrajectories(
    #     out_dirs,
    #     dir_name='output_',
    #     analysis_dir=str(Path(analysis_dir, "traj")),
    #     nproc=50
    # )
    #
    # test_dir = Path(Path.home(), "mtorc2/exp_10/94/analysis/test")
    # gsms_A_file = Path(test_dir, "test.csv")
    #
    # rmf_A_file = Path(test_dir, "test.rmf3")
    # scores_A_file = Path(test_dir, "test.txt")

    # HA = AT.get_models_to_extract(str(gsms_A_file))
    #
    # AT.extract_models_to_single_rmf(
    #     HA,
    #     str(rmf_A_file),
    #     str(output_dir),
    #     str(scores_A_file)
    # )

    output_dir = Path(Path.home(), "mtorc2/exp_10/94/sample/output.114524")
    analysis_dir = Path(Path.home(), "mtorc2/exp_10/94/analysis/176064")

    gsms_dir = Path(analysis_dir, "sampcon_0/gsms/5000")
    A_rmf3 = Path(gsms_dir, "A.rmf3")

    sampcon_dir = Path(analysis_dir, "sampcon_0/213807")
    for frame in [1,7,8,10,11,12,16,17,21,23,25,26,27,30,32]:
        get_rmsd(
            rmf_1=Path(sampcon_dir, "cluster.0/cluster_center_model.rmf3"),
            frame_1=0,
            rmf_2=Path(gsms_dir, "A.rmf3"),
            frame_2=frame
        )
