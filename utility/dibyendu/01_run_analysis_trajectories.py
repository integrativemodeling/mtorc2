# Analysis using PMI analysis
# (1) initalize the analysis class - AnalysisTrajectories (run_analysis_trajectories.py)
# (2) Add restraints you want to be analyzed
# (3) Read stat files to obtain the scores, nuisance parameters, and info about the rmf files


import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# It is assumed that reader knows about https://github.com/salilab/PMI_analysis/tree/main/pyext/src
#sys.path.append('/XX/XX/XX/imp_lib/PMI_analysis/pyext/src/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 16
top_dir = "./"
analys_dir = "./analysis"

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run' # if you have named your output folders from conformational sampling differently, change this line
out_dirs = glob.glob(top_dir+ dir_head+ '*/output/')
print(out_dirs)



################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'NHSF':30.0}     # this could be system dependent 

# 1. Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc,
                          nskip=20)
AT.ambiguous_XLs_restraint= True
AT.Multiple_XLs_restraints= True

# 2. Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
                             ambiguous_XLs_restraint = True,
                             Multiple_XLs_restraints = False,
                             get_nuisances = True)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()

# 3. Read stat files
AT.read_stat_files()
AT.write_models_info()
