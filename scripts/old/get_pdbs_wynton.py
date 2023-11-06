import os
from pathlib import Path
import time


def get_pdb(
        path_home,
        exp_num,
        traj_num,
        job_num
):
    pdb_dir = Path("/wynton/home/sali/mhancock/mtorc2/single_traj/exp_{}/{}/sample/output.{}/output_0/pdbs".format(exp_num, traj_num, job_num))
    pdb_file = Path(pdb_dir, "model.0.pdb")

    dest_dir = Path(path_home, "Downloads")
    dest_file = Path(dest_dir, "{}.{}.pdb".format(traj_num, job_num))

    scp_command = "scp mhancock@sali-log1.wynton.ucsf.edu:{} {}".format(pdb_file, dest_file)
    # print(scp_command)
    os.system(scp_command)


def get_output(
        path_home,
        exp_num,
        traj_num,
        job_num
):
    output_dir = Path("/wynton/home/sali/mhancock/mtorc2/exp_{}/{}/sample/output.{}".format(exp_num, traj_num, job_num))

    dest_dir = Path(path_home, "mtorc2/exp_{}/{}/sample".format(exp_num, traj_num))

    scp_command = "rsync -r mhancock@sali-log1.wynton.ucsf.edu:{} {}".format(output_dir, dest_dir)
    # print(scp_command)
    os.system(scp_command)


if __name__ == "__main__":
    exp_infos = list()
    exp_infos.append((8, 65, 90001, ["57", "A", "akt"]))
    # exp_infos.append((8, 66, 90004, ["57", "A", "akt"]))
    # exp_infos.append((8, 68, 81259, ["57", "A", "B", "akt"]))
    # exp_infos.append((8, 69, 81261, ["57", "A", "B", "akt"]))
    # exp_infos.append((8, 71, 81269, ["57", "A", "akt"]))

    for exp_info in exp_infos:
        t0 = time.time()
        get_output(
            path_home=Path.home(),
            exp_num=exp_info[0],
            traj_num=exp_info[1],
            job_num=exp_info[2]
        )
        print("{} took {}s".format(exp_infos[0][1], time.time()-t0))

        # get_pdb(
        #     path_home=Path.home(),
        #     exp_num=exp_info[0],
        #     traj_num=exp_info[1],
        #     job_num=exp_info[2]
        # )