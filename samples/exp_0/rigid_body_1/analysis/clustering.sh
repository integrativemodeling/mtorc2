#!/bin/bash
#$ -q gpu.q
#$ -N MTORC_cluster
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l mem_free=30G
#$ -l scratch=25G
#$ -l compute_cap=80,gpu_mem=23G

module load Sali
module load cuda/9.1
module load python3/pyrmsd/4.1.gita558b8a
module load imp-fast/2.14.0

top_dir="/wynton/home/sali/mhancock/mtorc2/rigid_body_1"
analys_dir="/wynton/home/sali/mhancock/mtorc2/rigid_body_1/data/analysis/cluster1"
name=MTORC
export top_dir
export analys_dir
export name

echo "START"
python3 /wynton/home/sali/mhancock/imp-sampcon/pyext/src/exhaust.py --sysname $name --path $analys_dir --mode cuda --cores 4 --align --density density.txt --gridsize 5.0 --gnuplot --scoreA A_models_score1.txt --scoreB B_models_score1.txt --rmfA A_models_clust1.rmf3 --rmfB B_models_clust1.rmf3
echo "END"
