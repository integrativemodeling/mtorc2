#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=96:00:00
#$ -l mem_free=5G
#$ -l scratch=200G
#$ -N sample_115
#$ -t 1-5
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

#cd "./test"
#LOC_OUT_DIR="./output_0"

N_FRAMES=20000
XLS=all_dss:1:DSS,all_edc:1:EDC
MAP=J238_3/J238_3_dsfact2_cutoff.01_ng200
ACTIVE_SITE=0
EM_COMPS=CRIM,RBD,MSIN1PH,KINASEN,KINASEC

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_114.py .

# Copy all of the subunit gmms (.mrc).
cp "$HOME"/mtorc2/data/em/comps/comps_114_5/* .

mpirun -np 8 python3 sample_114.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$MAP" "$ACTIVE_SITE" "$EM_COMPS"

# Move the local output into the global output directory.
mv "$LOC_OUT_DIR" "$GLOB_OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT

#cd ..