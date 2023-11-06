#!/bin/bash

#$ -S /usr/bin/python3
#$ -cwd
#$ -l h_rt=300:00:00
#$ -l mem_free=24G
#$ -l scratch=50G
#$ -N RICTOR_cpu
#$ -l hostname='qb3-id*'
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
# Adapted from alphafold/docker/run_alphafold.py script.
# Original version runs AlphaFold using a docker image.
# This adapted version uses a singularity image with defaults
# set for the UCSF Wynton cluster.


python run_alphafold21_cpu.py --fasta_paths=/wynton/home/sali/mhancock/mtorc2/data/fasta/RICTOR.fasta --model_preset=multimer --is_prokaryote_list=false
