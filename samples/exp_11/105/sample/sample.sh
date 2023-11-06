#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_105
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

N_FRAMES=500
XLS=all_domain_dss:1:DSS,all_domain_edc:1:EDC
CONN=15
MAP=J238_15_dsfact2_cutoff.1_ng200
EM_COMPS=CRIM,KINASE
RES_PER_COMP=5
REMOVED_REGIONS=MLST8:325-326
KAPPA=2
RADIUS=2

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_100.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_"$RES_PER_COMP"/* .

mpirun -np 8 python3 sample_100.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM_COMPS" "$RES_PER_COMP" "$MAP" "$REMOVED_REGIONS" "$KAPPA" "$RADIUS"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT