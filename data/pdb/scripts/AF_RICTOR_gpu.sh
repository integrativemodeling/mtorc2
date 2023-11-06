#!/bin/bash

#$ -S /usr/bin/python3
#$ -q gpu.q
#$ -N RICTOR_gpu
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l mem_free=60G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=46G
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
# Adapted from alphafold/docker/run_alphafold.py script.
# Original version runs AlphaFold using a docker image.
# This adapted version uses a singularity image with defaults
# set for the UCSF Wynton cluster.


python run_alphafold21_cpu.py --fasta_paths=/wynton/home/sali/mhancock/mtorc2/data/RICTOR.fasta --model_preset=multimer --is_prokaryote_list=false
