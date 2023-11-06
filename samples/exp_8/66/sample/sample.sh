#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=72:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N MTORC_66
#$ -t 1-10
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_65.py "$TMPDIR"
cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
GLOB_OUT_DIR="$SGE_CWD_PATH"/output."$JOB_ID"
mkdir "$GLOB_OUT_DIR"

N_FRAMES=100000
XLS=57,A,akt
CONN=20
EM=1
FULL_FLEX=0

echo N_FRAMES: "$N_FRAMES"
echo XLS: "$XLS"
echo CONN: "$CONN"
echo EM: "$EM"
echo FULL_FLEX: "$FULL_FLEX"

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_1/* .

mpirun -np 8 python3 sample_65.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM" "$FULL_FLEX"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT