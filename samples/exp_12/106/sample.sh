#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_106
#$ -t 1-1
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

N_FRAMES=2000
MAP=J238_3/J238_3_dsfact2_cutoff.01_ng200
EM_COMPS=KINASE
RES_PER_COMP=1

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_106.py .

# Copy all of the subunit gmms (.mrc).
#cp "$HOME"/mtorc2/data/em/components_"$RES_PER_COMP"/* .

mpirun -np 8 python3 sample_106.py "$LOC_OUT_DIR" "$N_FRAMES" "$EM_COMPS" "$MAP" "$RES_PER_COMP"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT