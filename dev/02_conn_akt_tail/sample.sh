#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N conn
#$ -t 1-1
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_217

N_FRAMES=1000

mpirun -np 1 python3 test_conn.py "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT