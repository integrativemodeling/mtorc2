import pandas as pd
from pathlib import Path
import IMP
import IMP.rmf
import RMF


def add_secondary_structure(
        ss_file,
        orig_file,
        merged_file
):
    # Reading data from file1
    with open(ss_file) as fp:
        data = fp.read()

    # Reading data from file2
    with open(orig_file) as fp:
        data2 = fp.read()

    # Merging 2 files
    # To add the data of file2
    # from next line
    data += "\n"
    data += data2

    with open(merged_file, 'w') as fp:
        fp.write(data)


def correct_chains(
        orig_file,
        new_file,
        missing_chains
):
    max_res_id_dict = {"A": 2549, "B": 1707, "C": 326, "D": 522, "E": 480, "F": 2549, "G": 1707, "H": 326, "I": 522, "J": 480}

    f_in = open(orig_file, "r")
    f_out = open(new_file, "w")
    f_in.readline()
    x = f_in.readline()
    chain = "A"
    while "ENDMDL" not in x:
        if "UNK" in x:
            x = x[:21] + chain + x[22:]
            x = x[:17] + "BEA" + x[20:]

        # print(x)
        f_out.write(x)

        res_id = int(x[22:26])
        if res_id == max_res_id_dict[chain]:
            if ord(chain) < 70:
                chain = chr(ord(chain) + 5)
            else:
                chain = chr(ord(chain) - 4)

            if missing_chains and chain in missing_chains:
                chain = chr(ord(chain) + 1)

        x = f_in.readline()

    f_out.write("ENDMDL")
    f_in.close()
    f_out.close()

    # Clean up the original pdb file.
    # orig_file.unlink()


def rmf_to_pdb(
        rmf_file,
        frame,
        pdb_file,
        # cluster_dir,
        missing_chains=None
):
    # rmf_file = Path(cluster_dir, "cluster_center_model.rmf3")
    print(rmf_file)

    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, frame)

    # pdb_file = Path(cluster_dir, "cluster_center_model.raw.pdb")
    # fixed_file = Path(cluster_dir, "cluster_center_model.pdb")
    tmp1_file = Path(str(pdb_file) + ".tmp1")
    tmp2_file = Path(str(pdb_file) + ".tmp2")
    IMP.atom.write_pdb(h, str(tmp1_file))
    print(pdb_file)

    # correct_chains(
    #     orig_file=tmp1_file,
    #     new_file=tmp2_file,
    #     missing_chains=missing_chains
    # )
    #
    # add_secondary_structure(
    #     ss_file=Path(Path.home(), "mtorc2/data/pdb/secondary_structure_exp.pdb"),
    #     orig_file=tmp2_file,
    #     merged_file=pdb_file
    # )

    # tmp1_file.unlink()
    # tmp2_file.unlink()


if __name__ == "__main__":
    home_dir = Path("/wynton/home/sali/mhancock")
    mtorc2_dir = Path(home_dir, "mtorc2")
    analysis_dir = Path(mtorc2_dir, "exp_10/89/analysis")
    sampcon_dir = Path(analysis_dir, "2462/sampcon_1/9510")
    cluster_dirs = [cluster_dir for cluster_dir in sampcon_dir.glob("*") if cluster_dir.is_dir()]
    for cluster_dir in cluster_dirs:
        print(cluster_dir)
        rmf_to_pdb(
            rmf_file=Path(cluster_dir, "cluster_center_model.rmf3"),
            frame=0,
            pdb_file=Path(cluster_dir, "cluster_center_model.pdb")
        )

    # rmf_file = Path(analysis_dir, "218765/sampcon_0/gsms/1000/A.rmf3")
    # for i in range(25):
    #     pdb_file = Path(analysis_dir, "218765/sampcon_0/gsms/1000/A_{}.pdb".format(i))
    #     rmf_to_pdb(
    #         rmf_file=rmf_file,
    #         frame=i,
    #         pdb_file=pdb_file
    #     )

    # params = list()
    # params.append(("exp_9/83", "817731", "sampcon_0.0", 0))
    # params.append(("exp_9/83", "817731", "sampcon_0.0", 1))
    #
    # for param in params:
    #     print(param)
    #     rmf_to_pdb(
    #         exp_dir=param[0],
    #         job_id=param[1],
    #         sampcon=param[2],
    #         c=param[3]
    #     )

