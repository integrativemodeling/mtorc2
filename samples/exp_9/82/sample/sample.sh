#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_82
#$ -t 1
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_217


cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
GLOB_OUT_DIR="$SGE_CWD_PATH"/output."$JOB_ID"
mkdir "$GLOB_OUT_DIR"

N_FRAMES=2500
XLS=81,A,B
CONN=25
EM=0
EM_COMPS=
AKT_DR=1
AKT_PH=1
SHUFFLE=1
PARTIAL_FLEX=1
MOVE_AWARE=0

echo N_FRAMES:     "$N_FRAMES"
echo XLS:          "$XLS"
echo CONN:         "$CONN"
echo EM:           "$EM"
echo EM_COMPS:     "$EM_COMPS"
echo AKT_DR:       "$AKT_DR"
echo AKT_PH:       "$AKT_PH"
echo SHUFFLE:      "$SHUFFLE"
echo PARTIAL_FLEX: "$PARTIAL_FLEX"
echo MOVE_AWARE:   "$MOVE_AWARE"

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_71.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_2/* .

mpirun -np 1 python3 sample_71.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM" "$EM_COMPS" "$AKT_DR" "$AKT_PH" "$SHUFFLE" "$PARTIAL_FLEX" "$MOVE_AWARE"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT