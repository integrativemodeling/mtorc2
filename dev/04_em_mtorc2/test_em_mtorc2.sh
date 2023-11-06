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
MAP="J238_3_mtorc2/J238_3_mtorc2_dsfact2_cutoff.05_ng100"
mpirun -np 8 python3 test_em_mtorc2.py "$N_FRAMES" "$MAP"