import pandas as pd
from pathlib import Path
import IMP
import IMP.rmf
import RMF
import math
import shutil


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


def update_chain(
        x_1,
        x_2
):
    norm = math.sqrt((x_1[0] - x_2[0])**2 + (x_1[1] - x_2[1])**2 + (x_1[2] - x_2[2])**2)

    # Arbitrary cutoff
    if norm > 25:
        return True
    else:
        return False


def next_chain(
    chain
):
    if ord(chain) < 70:
        new_chain = chr(ord(chain) + 5)
    else:
        new_chain = chr(ord(chain) - 4)

    return new_chain


def correct_chains(
        orig_file,
        new_file,
        missing_chains
):
    max_res_id_dict = {"A": 2549, "B": 1707, "C": 326, "D": 522, "E": 480}

    f_in = open(orig_file, "r")
    f_out = open(new_file, "w")
    # Burn first line
    f_in.readline()
    x = f_in.readline()
    chain = "A"
    while "ENDMDL" not in x:
        if "UNK" in x:
            x = x[:21] + chain + x[22:]
            x = x[:17] + "BEA" + x[20:]

        # print(x)
        f_out.write(x)

        # res_id = int(x[22:26])
        # if res_id == max_res_id_dict[chain]:
        #     if ord(chain) < 70:
        #         chain = chr(ord(chain) + 5)
        #     else:
        #         chain = chr(ord(chain) - 4)
        #
        #     if missing_chains and chain in missing_chains:
        #         chain = chr(ord(chain) + 1)

        res_id_prev = int(x[22:26])
        # c_prev = float(x[31:38]), float(x[39:46]), float(x[47:54])
        x = f_in.readline()

        if "ENDMDL" not in x:
            # c = float(x[31:38]), float(x[39:46]), float(x[47:54])
            res_id = int(x[22:26])

            if (res_id_prev > 10) and (res_id < 10):
                chain = next_chain(chain)

    f_out.write("ENDMDL")
    f_in.close()
    f_out.close()

    # Clean up the original pdb file.
    # orig_file.unlink()


def rmf_to_pdb(
        rmf_file,
        frame,
        pdb_file,
        corr_chains=True,
        add_ss=True,
        missing_chains=None
):
    print(rmf_file)

    fh = RMF.open_rmf_file_read_only(str(rmf_file))
    m = IMP.Model()
    h = IMP.rmf.create_hierarchies(fh, m)[0]
    IMP.rmf.load_frame(fh, frame)

    tmp_file = Path(str(pdb_file) + ".tmp")
    IMP.atom.write_pdb(h, str(pdb_file))
    print(pdb_file)

    if corr_chains:
        shutil.copy(pdb_file, tmp_file)
        correct_chains(
            orig_file=tmp_file,
            new_file=pdb_file,
            missing_chains=missing_chains
        )
        tmp_file.unlink()

    if add_ss:
        shutil.copy(pdb_file, tmp_file)
        add_secondary_structure(
            ss_file=Path(Path.home(), "mtorc2/data/pdb/secondary_structure_exp.pdb"),
            orig_file=tmp_file,
            merged_file=pdb_file
        )
        tmp_file.unlink()


if __name__ == "__main__":
    # cluster_dir = Path("/wynton/group/sali/mhancock/mtorc2/samples/exp_14/126/analysis/0/sampcon_2/3/cluster.0")
    cluster_dir = Path(Path.home(), "mtorc2/manuscript/submission_2/models/dev_09_2_-1_2/cluster.0")
    print(cluster_dir)
    rmf_to_pdb(
        rmf_file=Path(cluster_dir, "cluster_center_model.rmf3"),
        frame=0,
        pdb_file=Path(cluster_dir, "cluster_center_model.pdb"),
        add_ss=True,
        corr_chains=True
    )
