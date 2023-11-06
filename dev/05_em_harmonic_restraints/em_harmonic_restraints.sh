#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=8:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N test_em
#$ -t 1-1
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_217

N_FRAMES=1000
mpirun -np 1 python3 em_harmonic_restraints.py "$N_FRAMES"