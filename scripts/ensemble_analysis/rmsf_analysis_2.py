from pathlib import Path
import numpy as np

if __name__ == "__main__":
    job_dir = Path(Path.home(), "mtorc2/manuscript/submission_2/models/126_0_2_3/")
    rmsf_file = Path(job_dir, "cluster.0.rmsf.npy")
    avg_file = Path(job_dir, "cluster.0.avg.npy")

    mat_dir = Path("/wynton/group/sali/mhancock/mtorc2/tmp")
    mat_files = [mat_file for mat_file in mat_dir.glob("dist_*")]
    print(len(mat_files))
    dist_mat = np.load(mat_files[0])
    avg_mat = np.zeros(dist_mat.shape)
    # avg_mat.shape

    for mat_file in mat_files:
        print(mat_file)
        dist_mat = np.load(mat_file)
        avg_mat = avg_mat + dist_mat

        # break

    avg_mat = avg_mat / len(mat_files)
    np.save(avg_file, avg_mat)

    rmsf_mat = np.zeros(dist_mat.shape)
    # rmsf_mat.shape

    for mat_file in mat_files:
        print(mat_file)
        dist_mat = np.load(mat_file)
        rmsf_mat = rmsf_mat + (dist_mat - avg_mat)**2

        # break

    rmsf_mat = rmsf_mat / len(mat_files)
    rmsf_mat = rmsf_mat**(1/2)

    np.save(rmsf_file, rmsf_mat)

