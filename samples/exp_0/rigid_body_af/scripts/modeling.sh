#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=72:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N MTORC
#$ -pe smp 40
#$ -t 1-8
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
conda activate imp
module load Sali
echo $NSLOTS
mpirun -np $NSLOTS python modeling.py $SGE_TASK_ID 10000
