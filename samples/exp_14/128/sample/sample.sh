#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=72:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_128
#$ -t 1-10
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

EM_COMPS=CRIM,KINASEN,KINASEC
RES_PER_COMP=10
FLEX=1
XLS=dss:1:DSS,edc:1:EDC
MAP=J252/J252_4A/dsfact2_cutoff.05_ng100.txt
EM_W=10
SHUFFLE=1
N_FRAMES=30000

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_122.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/122/122_$RES_PER_COMP/* .

mpirun -np 8 python3 sample_122.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT