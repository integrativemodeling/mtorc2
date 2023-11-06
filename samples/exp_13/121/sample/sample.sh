#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=48:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_121
#$ -t 1-5
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_217

cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
#GLOB_OUT_DIR="$SGE_CWD_PATH"/output."$JOB_ID"
PATH_TAIL="${SGE_CWD_PATH#/*/*/*/*/}"
GLOB_OUT_DIR="/wynton/group/sali/mhancock/$PATH_TAIL/output.$JOB_ID"
mkdir -p "$GLOB_OUT_DIR"

#cd "./test"
#LOC_OUT_DIR="./output_0"

RES_PER_COMP=5
AKT_TAIL_RES=10
EM_COMPS=CRIM,KINASEN,KINASEC
FULL_ARCH_EV=0
XLS=all_dss:1:DSS,all_edc:1:EDC
MAP=J238_3/J238_3_dsfact2_cutoff.01_ng100
EM_W=5
ACTIVE_SITE=0
N_FRAMES=10000

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_114.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/em_arch_$RES_PER_COMP/* .

mpirun -np 8 python3 sample_114.py "$LOC_OUT_DIR" "$RES_PER_COMP" "$AKT_TAIL_RES" "$EM_COMPS" "$FULL_ARCH_EV" "$XLS" "$MAP" "$EM_W"  "$ACTIVE_SITE" "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT

#cd ..