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


module load CBI conda-stage
conda activate imp_217
module load Sali python3/pyrmsd/4.1.gita558b8a
module load cuda/7

SAMPLE_NUM=935122
ANALYSIS_NUM=2
ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/analysis/$ANALYSIS_NUM"
SAMPLE_DIR="/wynton/group/sali/mhancock/mtorc2/dev/09_conn_akt_tail_2/sample/output.$SAMPLE_NUM"

# Run Sampcon.
JOB_ID=3
CLUSTER=-1
N_STRUCT=5000
FILTER="R"
GSMS_FOLDER="5000R"
COPY=
THRESH=

python ~/mtorc2/dev/09_conn_akt_tail_2/scripts/tail_sampcon.py "$ANALYS_DIR" "$SAMPLE_DIR" "$JOB_ID" "$CLUSTER" "$N_STRUCT" "$FILTER" "$GSMS_FOLDER" "$COPY" "$THRESH"
