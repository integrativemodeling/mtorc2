#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=72:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_92
#$ -t 1-50
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

N_FRAMES=50000
XLS=all.domain
CONN=25
MAP=J238_15_dsfact2_cutoff.1_ng100
EM_COMPS=CRIM,KINASE,MSIN1FLEX,KINASEFLEX
REMOVED_REGIONS=MLST8:325-326
AKT1PH=1

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_86.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_10/* .

mpirun -np 8 python3 sample_86.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM_COMPS" "$MAP" "$REMOVED_REGIONS" "$AKT1PH"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT