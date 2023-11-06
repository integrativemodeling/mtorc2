#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N mTORC_71_test
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp

# Copy the modeling script into the job scratch dir.
LOC_OUT_DIR="./output"
N_FRAMES=100
XLS=all.domain
CONN=25
EM=1

# Copy all of the subunit gmms (.mrc).
cp $HOME/mtorc2/scripts/sample_71.py .
cp "$HOME"/mtorc2/data/em/components_2/* .

python3 sample_71.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT