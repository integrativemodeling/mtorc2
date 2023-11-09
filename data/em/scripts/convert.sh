#!/bin/bash


export map_dir_name=1031
export map=1031_10.mrc
export NG=200
export dsfact=2
export cutoff=.05

export map_dir="$HOME/mtorc2/data/em/maps/${map_dir_name}"
export map_file="${map_dir}/${map}"
export param_str="dsfact${dsfact}_cutoff${cutoff}_ng${NG}"
export gmm_gmm_file="${map_dir}/${param_str}.gmm"
export gmm_txt_file="${map_dir}/${param_str}.txt"
export gmm_mrc_file="${map_dir}/${param_str}.mrc"

~/gmconvert/gmconvert V2G -imap "${map_file}" -ng ${NG} -cutoff "$cutoff" -fdsmap ${dsfact} -max_memory 100000 -emalg DG -ogmm "${gmm_gmm_file}"

python write_IMP_file_2.py "${gmm_gmm_file}"
python txt_to_mrc.py "${gmm_txt_file}" "${map_file}"