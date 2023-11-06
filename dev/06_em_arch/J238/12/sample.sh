#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N em_arch
#$ -t 1-3
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

OUT_DIR="$LOC_OUT_DIR"
RES_PER_COMP=5
MAP="J238_3/J238_3_dsfact2_cutoff.05_ng100"
SCORE_SS=0
N_FRAMES=25000

cp "$HOME"/mtorc2/data/em/comps/em_arch_$RES_PER_COMP/* .
cp "$HOME"/mtorc2/dev/06_em_arch/em_arch_10.py .

mpirun -np 8 python3 em_arch_10.py "$OUT_DIR" "$RES_PER_COMP" "$MAP" "$SCORE_SS" "$N_FRAMES"

mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT