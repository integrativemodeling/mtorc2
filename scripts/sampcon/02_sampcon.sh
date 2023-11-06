#!/bin/bash
#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/mtorc2/samples/exp_14/129/analysis
#$ -N smpcn_129
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -q gpu.q
#$ -l hostname='qb3-id*'


#module load CBI conda-stage
#conda activate imp_217
#module load Sali python3/pyrmsd/4.1.gita558b8a
#module load cuda/7

EXP_NUM=14
TRAJ_NUM=129
SAMPLE_NUM=347466
ANALYSIS_NUM=1
ANALYS_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/analysis/$ANALYSIS_NUM"
SAMPLE_DIR="$HOME/mtorc2/samples/exp_$EXP_NUM/$TRAJ_NUM/sample/output.$SAMPLE_NUM"

# Run Sampcon.
JOB_ID=0
CLUSTER=3
N_STRUCT=10000
FILTER="R"
GSMS_FOLDER="10000R_edit"
RANGE=9
COPY=
THRESH=

python ~/mtorc2/scripts/sampcon/sampcon.py "$ANALYS_DIR" "$SAMPLE_DIR" "$JOB_ID" "$CLUSTER" "$N_STRUCT" "$FILTER" "$GSMS_FOLDER" "$RANGE" "$COPY" "$THRESH"
