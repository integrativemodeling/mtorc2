#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/samples/tmp/$JOB_ID.$TASK_ID.o
#$ -cwd
#$ -j y
#$ -l h_rt=18:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N mtor_130
#$ -t 1-100
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
# eval "$(conda shell.bash hook)"
# module load CBI conda-stage
# conda activate imp_219

module load Sali
module load imp
module load mpi
module load python3/pandas/0.25.3

EXP_ID=15
SAMPLE_ID=130

cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
# PATH_TAIL="${SGE_CWD_PATH#/*/*/*/*/}"
GLOB_OUT_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_ID/$SAMPLE_ID/$JOB_ID"
mkdir -p "$GLOB_OUT_DIR"

EM_COMPS=CRIM,KINASEN,KINASEC
RES_PER_COMP=10
FLEX=1
XLS=dss:1:DSS,edc:1:EDC
MAP=1031/dsfact2_cutoff.05_ng100.txt
EM_W=1
SHUFFLE=0
N_FRAMES=7500

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/130/10/* .

mpirun -np 8 python3 sample.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"
# python3 sample.py "$LOC_OUT_DIR" "$EM_COMPS" "$RES_PER_COMP"  "$FLEX" "$XLS" "$MAP" "$EM_W" "$SHUFFLE" "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT