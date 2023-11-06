from pathlib import Path
import sys
import time
import os
import shutil

sys.path.append(str(Path(Path.home(), "mtorc2/scripts/sampcon")))
import sampcon


def tail_sampcon(
        analysis_dir,
        gsms_dir,
        mode,
        thresh
):
    print("**Running sampcon***")
    t0 = time.time()

    rmf_A_file = Path(gsms_dir, "A.rmf3")
    rmf_B_file = Path(gsms_dir, "B.rmf3")
    scores_A_file = Path(gsms_dir, "A.txt")
    scores_B_file = Path(gsms_dir, "B.txt")

    custom_ranges_dir = Path(Path.home(), "mtorc2/dev/09_conn_akt_tail_2/data")
    selection_file = Path(custom_ranges_dir, "custom_range.txt")
    density_file = selection_file

    if mode == "cpu":
        mode_arg = "cpu_omp"
    elif mode == "gpu":
        mode_arg = "cuda"

    thresh_arg = ""
    if thresh:
        thresh_arg = "--skip -ct {}".format(thresh)

    cmd = "python3 {}/imp-sampcon/pyext/src/exhaust.py --sysname mtorc2 --path {} --rmfA {} --rmfB {} --scoreA {} --scoreB {} --selection {} --density {} --gnuplot --align -m {} {} -c 4 -g 2.0 -gp".format(
        Path.home(), analysis_dir, rmf_A_file, rmf_B_file, scores_A_file, scores_B_file, selection_file, density_file, mode_arg, thresh_arg)

    print(cmd)
    os.system(cmd)

    print("Finished in {}s".format(time.time()-t0))


if __name__ == "__main__":
    analysis_dir = Path(sys.argv[1])
    sample_dir = Path(sys.argv[2])
    job_id = os.environ["JOB_ID"] if sys.argv[3] == "" else sys.argv[3]
    cluster = int(sys.argv[4])
    N_struct = int(sys.argv[5])
    filter = sys.argv[6]
    gsms_folder = sys.argv[7]
    copy = None if sys.argv[8] == "" else sys.argv[8]
    thresh = None if sys.argv[9] == "" else int(sys.argv[9])

    mode = "gpu"

    print("analysis_dir: {}".format(analysis_dir))
    print("sample_dir:   {}".format(sample_dir))
    print("job_id:       {}".format(job_id))
    print("cluster:      {}".format(cluster))
    print("N_struct:     {}".format(N_struct))
    print("filter:       {}".format(filter))
    print("gsms_folder:  {}".format(gsms_folder))
    print("copy:         {}".format(copy))
    print("thresh:       {}".format(thresh))
    print("mode:         {}".format(mode))

    sampcon_dir = Path(analysis_dir, "sampcon_{}".format(cluster))
    if not sampcon_dir.exists():
        sampcon_dir.mkdir()
        Path(sampcon_dir, "gsms").mkdir()

    job_dir = Path(sampcon_dir, job_id)
    job_dir.mkdir()

    gsms_dir = Path(sampcon_dir, "gsms/{}".format(gsms_folder))

    if not gsms_dir.exists():
        gsms_dir.mkdir()
        sampcon.extract_models(
            output_dir=sample_dir,
            analysis_dir=analysis_dir,
            filter=filter,
            gsms_dir=gsms_dir,
            N=N_struct,
            cluster=cluster
        )
    else:
        print("gsms dir found")

    # # Need to copy over 3 files from an older sampcon run. The pre-computed distance matrix is used in clustering. Sampling precision statistics and Chi square grid statistics are for constructing plots.
    if copy:
        copy_files = list()
        copy_files.append("Distances_Matrix.data.npy")
        copy_files.append("mtorc2.Sampling_Precision_Stats.txt")
        copy_files.append("mtorc2.ChiSquare_Grid_Stats.txt")

        for file_name in copy_files:
            shutil.copy(Path(sampcon_dir, copy, file_name), Path(job_dir, file_name))

    os.chdir(job_dir)
    tail_sampcon(
        analysis_dir=analysis_dir,
        gsms_dir=gsms_dir,
        mode=mode,
        thresh=thresh
    )