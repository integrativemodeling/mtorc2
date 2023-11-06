#$ -S /bin/bash
#$ -cwd
#$ -o ./sample
#$ -j y
#$ -l h_rt=48:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N MTORC_60
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

cp modeling.py "$TMPDIR"
cd "$TMPDIR"
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp
mpirun -np 8 python3 modeling.py
[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT