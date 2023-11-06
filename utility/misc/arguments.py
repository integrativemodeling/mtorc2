import argparse
from pathlib import Path


def sample_command_line_arguments():
    # Setup the command line arguments
    parser = argparse.ArgumentParser(description='Script runs replica exchange sampling of MTORC2 complex')
    parser.add_argument('--run', type=int,
                        help='An optional integer argument')
    parser.add_argument('--n_steps', type=int,
                        help='Number of steps to sample ')
    parser.add_argument('--job_dir', type=str,
                        help='Sample on hpc. Writes output files to scratch')
    args = parser.parse_args()

    n_steps = args.n_steps
    run_id = args.run
    job_dir = Path(args.job_dir)

    return n_steps, run_id, job_dir