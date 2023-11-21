module load Sali
module load imp
module load mpi
module load python3/pandas/0.25.3


LOC_OUT_DIR=/wynton/home/sali/mhancock/mtorc2/samples/test/output_0
rm -r "$LOC_OUT_DIR"
mkdir "$LOC_OUT_DIR"

cd /wynton/home/sali/mhancock/mtorc2/samples/test

PARAM_FILE=/wynton/home/sali/mhancock/mtorc2/data/params/130.csv
EM_COMPS=CRIM,KINASEN,KINASEC
RES_PER_COMP=10
SYM=1
FLEX=1
XLS=dss:1:DSS,edc:1:EDC
MAP=/wynton/home/sali/mhancock/mtorc2/data/em/maps/1031/1031_10_crop_dsfact2_cutoff.05_ng100.txt
EM_W=1
SHUFFLE=0
N_FRAMES=100

# Copy the modeling script into the job scratch dir.
cp $HOME/mtorc2/scripts/sample/sample.py .

python3 sample.py --output_dir "$LOC_OUT_DIR" --param_file "$PARAM_FILE" --em_comps "$EM_COMPS" --res_per_comp "$RES_PER_COMP" --sym "$SYM" --flex "$FLEX" --xls "$XLS" --map "$MAP" --em_w "$EM_W" --shuffle "$SHUFFLE" --n_frames "$N_FRAMES"

cd /wynton/home/sali/mhancock/mtorc2/scripts/sample