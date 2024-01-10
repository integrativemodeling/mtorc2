#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/samples/tmp/test.$JOB_ID.$TASK_ID.o
#$ -cwd
#$ -j y
#$ -l h_rt=00:5:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N test_mpi
#$ -t 1-1
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'


eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219

mpirun -np 8 python $HOME/mtorc2/dev/10_mpi/test.py