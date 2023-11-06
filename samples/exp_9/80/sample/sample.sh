#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_80
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

N_FRAMES=25000
XLS=all.domain
CONN=25
EM=1
EM_COMPS=MTOR,RICTOR,MSIN1,MLST8,CRIM,KINASE
NO_COPY=
AKT_DR=0

echo N_FRAMES: "$N_FRAMES"
echo XLS: "$XLS"
echo CONN: "$CONN"
echo EM: "$EM"
echo EM_COMPS: "$EM_COMPS"
echo NO_COPY: "$NO_COPY"
echo AKT_DR: "$AKT_DR"

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_71.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/components_2/* .

mpirun -np 8 python3 sample_71.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM" "$EM_COMPS" "$NO_COPY" "$AKT_DR"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT