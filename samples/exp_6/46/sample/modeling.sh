#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=72:00:00
#$ -l mem_free=25G
#$ -l scratch=25G
#$ -N MTORC_46
#$ -pe smp 8
#$ -l hostname='qb3-id*'

echo "$TMPDIR"
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp
echo $NSLOTS
mpirun -np $NSLOTS python3 modeling.py "$TMPDIR"
[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT