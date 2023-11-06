from pathlib import Path
import random


def condense_outputs(
        output_dirs
):
    output_dict = {"run": list(), "tmp": list()}

    tmp_output_dir = Path(
        "/wynton/scratch/mhancock/output{}".format(random.random()))
    os.system("mkdir {}".format(tmp_output_dir))
    print(tmp_output_dir)
    run_num = 0
    for output_dir in output_dirs:
        run_dirs = [run_dir for run_dir in output_dir.glob("*") if
                    run_dir.is_dir()]
        for run_dir in run_dirs:
            tmp_run_dir = Path(tmp_output_dir, "output_{}".format(run_num))

            output_dict["run"].append(run_dir)
            output_dict["tmp"].append(tmp_run_dir)
            os.system("cp -r {} {}".format(run_dir, tmp_run_dir))

            run_num = run_num + 1

    output_df = pd.DataFrame(output_dict)
    return output_df, tmp_output_dir