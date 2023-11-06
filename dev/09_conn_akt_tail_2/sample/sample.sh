#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=16:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N tail
#$ -t 1-25
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218

cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
PATH_TAIL="${SGE_CWD_PATH#/*/*/*/*/}"
GLOB_OUT_DIR="/wynton/group/sali/mhancock/$PATH_TAIL/output.$JOB_ID"
mkdir -p "$GLOB_OUT_DIR"

PDB_FILE="$HOME/mtorc2/manuscript/submission_2/models/126_0_2_3/cluster.0/cluster_center_model.pdb"
XLS=dss:1:DSS,edc:1:EDC
PHOS_RES=450
SHUFFLE=1
N_FRAMES=25000

echo "LOC_OUT_DIR: $LOC_OUT_DIR"
echo "PDB_FILE: $PDB_FILE"
echo "XLS: $XLS"
echo "PHOS_RES: $PHOS_RES"
echo "SHUFFLE: $SHUFFLE"
echo "N_FRAMES: $N_FRAMES"

cp $HOME/mtorc2/dev/09_conn_akt_tail_2/scripts/model_akt_tail.py .

mpirun -np 8 python3 model_akt_tail.py "$LOC_OUT_DIR" "$PDB_FILE" "$XLS" "$PHOS_RES" "$SHUFFLE" "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT