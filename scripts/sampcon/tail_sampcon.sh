#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis
#$ -N tail_sampcon
#$ -j y
#$ -l h_rt=4:00:00
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

SAMPLE_NUM=34291228
ANALYSIS_NUM=1
ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis/$ANALYSIS_NUM"
SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/sample/$SAMPLE_NUM"

# Run Sampcon.
JOB_ID=1
SCORE_CLUSTER=-1
N_STRUCT=1000
FILTER="R"
GSMS_DIR_NAME="1000R"
RANGE_FILE="/wynton/home/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/data/custom_range.txt"
COPY=0
THRESH=50

# python ~/mtorc2/scripts/sampcon/sampcon.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --job_id "$JOB_ID" --score_cluster "$SCORE_CLUSTER" --n_struct "$N_STRUCT" --filter "$FILTER" --gsms_dir_name "$GSMS_DIR_NAME" --range_file "$RANGE_FILE"

python ~/mtorc2/scripts/sampcon/sampcon.py --analysis_dir "$ANALYS_DIR" --sample_dir "$SAMPLE_DIR" --job_id "$JOB_ID" --score_cluster "$SCORE_CLUSTER" --n_struct "$N_STRUCT" --filter "$FILTER" --gsms_dir_name "$GSMS_DIR_NAME" --range_file "$RANGE_FILE" --copy "$COPY" --thresh "$THRESH"