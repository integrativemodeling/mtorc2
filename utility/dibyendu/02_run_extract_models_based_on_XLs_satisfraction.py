#!/usr/bin/env python3

import IMP
import IMP.rmf
import RMF
import pandas as pd
import multiprocessing as mp
import glob
import sys
import shutil
import re
import math
import os

# concept of the codes are taken from https://github.com/salilab/PMI_analysis/blob/main/pyext/src/analysis_trajectories.py

def extract_models(run_id, XL_sat_sorted_dict, all_rmf3_top_dir, analysis_path, gsms_dir):
        '''
        Use rmf_slice to extract the GSMs
        '''
        num = re.search(r'\d+', run_id).group()
        csv_file = 'scores_info_' + num + '.csv'
        gsms_info =pd.read_csv(os.path.join(analysis_path, csv_file))
        
        all_frames_to_extract = XL_sat_sorted_dict[run_id]
        for row in gsms_info.itertuples():
            ID = run_id
            fr = int(row.MC_frame)
            if fr in all_frames_to_extract:
                all_frames_to_extract.remove(fr)
                print('MC_Frame: ', run_id, fr)
                fr_rmf = int(row.rmf_frame_index)
                rmf_file = row.rmf3_file
                traj_in = os.path.join(all_rmf3_top_dir, rmf_file)
                if row.half == 'A':
                    filename = 'h1'
                else:
                    filename = 'h2'

                file_out = os.path.join(
                    gsms_dir, filename+'_'+str(ID)+'_'+str(fr)+'.rmf3')

                os.system('rmf_slice -q ' + traj_in + ' ' + file_out
                          + ' --frame ' + str(fr_rmf))
            else:
                pass


def do_extract_models(nproc, XL_sat_sorted_dict, analysis_path, gsms_dir):

    all_runs_ids = list(XL_sat_sorted_dict.keys())
    total_num_runs = len(all_runs_ids)
    num_iter = math.ceil(total_num_runs/nproc)

    
    data = []
    for i in range(total_num_runs):
        # Setup a list of processes that we want to run
        data.append((all_runs_ids[i], XL_sat_sorted_dict, all_rmf3_top_dir, analysis_path, gsms_dir)) 
    
    with mp.Pool(processes=nproc) as pool:
        pool.starmap(extract_models, data)


nproc = 12
analysis_path = '/XX/XX/XX/analysis' # This path depends on your path that you defined in 01_run_analysis_trajectories.py
gsms_dir = os.path.join(analysis_path, 'Good_models') # name however you want your new directory should be
all_rmf3_top_dir = '/XX/XX/XX'
isExist = os.path.exists(gsms_dir)

# 1. Make a folder to store all the models with atleast 80% XLs satisfaction
if isExist:
    shutil.rmtree(gsms_dir)

os.makedirs(gsms_dir)


# 2. This information can be obtained from other_info_*.csv file; these files should be in your analysis folder that you should have generated 
# after using 01_run_analysis_trajectories.py file

all_other_csv = glob.glob('other_info_*.csv')

os.chdir(analysis_path)

run_and_mc_frame = {}
list_df = []
# 3. check each csv files; one file from one independent conformational sampling
for f_csv in all_other_csv:
    temp_csv_num = f_csv.split('.')[0].split('_')[-1]
    temp_run = 'run'+str(temp_csv_num)
    df = pd.read_csv(f_csv)
    run_and_mc_frame[temp_run] = df.loc[df['XLs_satif'] >0.80, 'MC_frame'].tolist()

sorted_dict = dict(sorted(run_and_mc_frame.items()))

print("searching started .... ")

# 4. extract all the rmf3 files in a GSM (good scoring model folder)
do_extract_models(nproc, sorted_dict, analysis_path, gsms_dir)

print('done')
