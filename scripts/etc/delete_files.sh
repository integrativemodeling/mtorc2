#!/bin/bash


EXP_DIR=/wynton/group/sali/mhancock/mtorc2/samples/exp_13
for JOB_NUM in {114..121}
do
    JOB_DIR="$EXP_DIR/$JOB_NUM"
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &

    # for OUT_DIR in "$JOB_DIR"/output_*; do
    #     # Extract the Z value from the subfolder name using string manipulation
    #     folder_name=$(basename "$OUT_DIR")  # Get the base folder name without the path
    #     z_value=${folder_name#output_}        # Removes "output_" from the beginning of the folder name

    #     # Check if Z is greater than or equal to 50 using an if statement
    #     if [ "$z_value" -ge 50 ]; then
    #         echo "$OUT_DIR"
    #         rm -r "$OUT_DIR"
    #     fi
    # done

done