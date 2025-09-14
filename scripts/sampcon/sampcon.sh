#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/analysis/tmp/$JOB_ID.o
#$ -N smpcn_136
#$ -j y
#$ -l h_rt=00:10:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -q gpu.q
#$ -l hostname='qb3-id*'


# module load CBI conda-stage
# conda activate imp_217
# module load Sali
# module load python3/pyrmsd/4.1.gita558b8a
# module load cuda/7

# export CUDA_VISIBLE_DEVICES=1

EXP_NUM=15
TRAJ_NUM=136
SAMPLE_NUM=966247
ANALYSIS_NUM=2
ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/analysis/$TRAJ_NUM/$ANALYSIS_NUM"
SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/$SAMPLE_NUM"

# Run Sampcon.
JOB_ID=3
SCORE_CLUSTER=-1
N_STRUCT=1000
FILTER="R"
GSMS_DIR_NAME="1000R"
RANGE_FILE="/wynton/home/sali/mhancock/mtorc2/data/ranges/custom_ranges_ph.txt"
SYMM_FILE="/wynton/home/sali/mhancock/mtorc2/data/ranges/symm_groups_ph.txt"
COPY=2
THRESH=50

# python ~/mtorc2/scripts/sampcon/sampcon.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --job_id "$JOB_ID" --score_cluster "$SCORE_CLUSTER" --n_struct "$N_STRUCT" --filter "$FILTER" --gsms_dir_name "$GSMS_DIR_NAME" --range_file "$RANGE_FILE" --symm_file "$SYMM_FILE"

python ~/mtorc2/scripts/sampcon/sampcon.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --job_id "$JOB_ID" --score_cluster "$SCORE_CLUSTER" --n_struct "$N_STRUCT" --filter "$FILTER" --gsms_dir_name "$GSMS_DIR_NAME" --range_file "$RANGE_FILE" --symm_file "$SYMM_FILE" --copy "$COPY" --thresh "$THRESH"

# python ~/mtorc2/scripts/sampcon/sampcon.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --job_id "$JOB_ID" --score_cluster "$SCORE_CLUSTER" --n_struct "$N_STRUCT" --filter "$FILTER" --gsms_dir_name "$GSMS_DIR_NAME" --range_file "$RANGE_FILE"
