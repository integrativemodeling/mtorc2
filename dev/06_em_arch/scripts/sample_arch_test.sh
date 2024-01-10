module load Sali
module load imp
module load mpi
module load python3/pandas/0.25.3


LOC_OUT_DIR="/wynton/group/sali/mhancock/mtorc2/dev/06_em_arch/test/output_0"
mkdir -p "$LOC_OUT_DIR"

MAP=1031/dsfact2_cutoff.05_ng100.txt
N_FRAMES=1000

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/dev/06_em_arch/scripts/sample_arch.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/130/10/* .

python3 sample_arch.py "$LOC_OUT_DIR" "$MAP" "$N_FRAMES"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT