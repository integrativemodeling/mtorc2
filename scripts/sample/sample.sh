#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/mtorc2/samples/tmp/$JOB_ID.$TASK_ID.o
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=5G
#$ -l scratch=1G
#$ -N mtor_137
#$ -t 1-50
#$ -pe mpi_onehost 8
#$ -l hostname='qb3-id*'

module load Sali
module load imp
module load mpi
module load python3/pandas/0.25.3

EXP_ID=15
SAMPLE_ID=137

cd "$TMPDIR"
LOC_OUT_DIR="$TMPDIR"/output_"$((SGE_TASK_ID-1))"
GLOB_OUT_DIR="/wynton/group/sali/mhancock/mtorc2/samples/exp_$EXP_ID/$SAMPLE_ID/$JOB_ID"
mkdir "$LOC_OUT_DIR"
mkdir -p "$GLOB_OUT_DIR"

PARAM_FILE=/wynton/home/sali/mhancock/mtorc2/data/params/130.csv
RES_PER_COMP=10
SYM=1
FLEX=1
XLS=dss:1:DSS,edc:1:EDC
MAP=/wynton/home/sali/mhancock/mtorc2/data/em/maps/1031/1031_5_dsfact2_cutoff.05_ng100.txt
EM_W=10
EV=1
CONN=1
SHUFFLE=0
N_FRAMES=10000

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample.py .

mpirun -np 8 python3 sample.py --output_dir "$LOC_OUT_DIR" --param_file "$PARAM_FILE" --res_per_comp "$RES_PER_COMP" --flex "$FLEX" --sym "$SYM" --xl "$XLS" --map "$MAP" --em_w "$EM_W" --ev "$EV" --conn "$CONN" --shuffle "$SHUFFLE" --n_frames "$N_FRAMES"


# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT