# Setup the conda environment.
#eval "$(conda shell.bash hook)"
#module load CBI conda-stage
#conda activate imp_217

mkdir ./test
cd ./test
LOC_OUT_DIR=./output_0

N_FRAMES=100
XLS=81,A,B
CONN=25
EM=0
EM_COMPS=
AKT_DR=0
AKT_PH=0
SHUFFLE=1
PARTIAL_FLEX=1
MOVE_AWARE=0

echo N_FRAMES:     "$N_FRAMES"
echo XLS:          "$XLS"
echo CONN:         "$CONN"
echo EM:           "$EM"
echo EM_COMPS:     "$EM_COMPS"
echo AKT_DR:       "$AKT_DR"
echo AKT_PH:       "$AKT_PH"
echo SHUFFLE:      "$SHUFFLE"
echo PARTIAL_FLEX: "$PARTIAL_FLEX"
echo MOVE_AWARE:   "$MOVE_AWARE"

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample_71.py .

mpirun -np 1 python3 sample_71.py "$LOC_OUT_DIR" "$N_FRAMES" "$XLS" "$CONN" "$EM" "$EM_COMPS" "$AKT_DR" "$AKT_PH" "$SHUFFLE" "$PARTIAL_FLEX" "$MOVE_AWARE"