#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=8:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N MTORC_75
#$ -t 1-10
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp


cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
GLOB_OUT_DIR="$SGE_CWD_PATH"/output."$JOB_ID"
mkdir "$GLOB_OUT_DIR"

N_FRAMES=1500
XLS=NA
CONN=25
EM=1
EM_COMPS=CRIM,KINASE
EM_TEST=1

echo N_FRAMES: "$N_FRAMES"
echo XLS: "$XLS"
echo CONN: "$CONN"
echo EM: "$EM"
echo EM_COMPS: "$EM_COMPS"
echo EM_TEST: "$EM_TEST"

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_71.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_2/* .

mpirun -np 1 python3 sample_71.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM" "$EM_COMPS" "$EM_TEST"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT