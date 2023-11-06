import os
from pathlib import Path


def wynton_setup():
    tmp_dir = os.getenv('TMPDIR', '1')
    orig_dir = os.getenv('SGE_CWD_PATH', '1')
    job_id = os.getenv('JOB_ID', '1')
    task_id = os.getenv('SGE_TASK_ID', '1')
    task_id_0 = int(task_id) - 1
    print("ENV VARIABLES")
    print("tmp_dir:  ", tmp_dir)
    print("orig_dir: ", orig_dir)
    print("job_id:   ", job_id)
    print("task_id:  ", task_id)

    # The local output directory (on /scratch).
    loc_out_dir = Path(tmp_dir, "output.{}".format(job_id), "output_{}".format(task_id_0))
    print("loc_out_dir: ", loc_out_dir)

    # Create the global output directory (on /wynton).
    glob_out_dir = Path(orig_dir, "sample/output.{}".format(job_id))
    cmd_1 = "mkdir {}".format(glob_out_dir)
    print(cmd_1)
    os.system(cmd_1)

    # Copy the .mrc files for the component GMMs to the local directory (on /scratch).
    density_dir = Path(Path.home(), "mtorc2/data/em/components_1")
    cmd_2 = "cp {}/* .".format(density_dir)
    print(cmd_2)
    os.system(cmd_2)

    return loc_out_dir, glob_out_dir


def wynton_cleanup(
        loc_out_dir,
        glob_out_dir
):
    # Move the local output directory to the global output directory.
    cmd_3 = "mv {} {}".format(loc_out_dir, glob_out_dir)
    print(cmd_3)
    os.system(cmd_3)