#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -l mem_free=10G
#$ -l scratch=10G
#$ -N MTORC
#$ -pe smp 10
#$ -t 1-100
#$ -l hostname='qb3-id*'

echo "$TMPDIR"
#mkdir "$TMPDIR"
eval "$(conda shell.bash hook)"
#module load Sali imp/nightly mpi
conda activate imp
echo $NSLOTS
mpirun -np $NSLOTS python3 ../scripts/01_modeling.py --run_id $SGE_TASK_ID --n_steps 100000 --base_dir "$TMPDIR"
scp -r "$TMPDIR"/mtorc2/rigid_body_1/data/output/* /wynton/home/sali/mhancock/mtorc2/rigid_body_1/data/output
[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
