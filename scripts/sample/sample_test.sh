module load Sali
module load imp
module load mpi
module load python3/pandas/0.25.3
module load python3/scikit/0.21.3


LOC_OUT_DIR=/wynton/home/sali/mhancock/mtorc2/samples/test/output_0
rm -r "$LOC_OUT_DIR"
mkdir "$LOC_OUT_DIR"

cd /wynton/home/sali/mhancock/mtorc2/samples/test
rm /wynton/home/sali/mhancock/mtorc2/samples/test/*

PARAM_FILE="$HOME/mtorc2/data/params/130.csv"
N_FRAMES=0

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample.py .

python3 sample.py  --output_dir "$LOC_OUT_DIR" --param_file "$PARAM_FILE" --flex --n_frames "$N_FRAMES"


cd /wynton/home/sali/mhancock/mtorc2/scripts/sample