#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/samples/tmp/$JOB_ID.$TASK_ID.o
#$ -cwd
#$ -j y
#$ -l h_rt=1:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N em
#$ -t 1-1
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


cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
# PATH_TAIL="${SGE_CWD_PATH#/*/*/*/*/}"
GLOB_OUT_DIR="/wynton/group/sali/mhancock/mtorc2/dev/06_em_arch/1031/$JOB_ID"
mkdir -p "$GLOB_OUT_DIR"

MAP=1031/dsfact2_cutoff.05_ng100.txt
N_FRAMES=1000

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/130/10/* .

mpirun -np 8 python3 sample.py "$LOC_OUT_DIR" "$EM_COMPS" "$N_FRAMES"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT