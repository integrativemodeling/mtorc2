import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('/home/rakesh/WORK/Projects/SPOTSComplex/IMPModeling/AnalysisModules/PMI_analysis/pyext/src')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################
nproc = 32
top_dir =  "/salilab/park3/rakesh/Projects/3-SPOTSComplex/3-IntModeling/ClosedFinal100"
analys_dir = top_dir+'/analys/'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'modelingRun'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/jobresults/')

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'DSSO':35.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs, ambiguous_XLs_restraint = True, Multiple_XLs_restraints = False)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
AT.set_analyze_EM_restraint()
AT.set_analyze_score_only_restraint('MembraneRestraint_Score', 'Mem_score')

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.get_psi_stats()

AT.hdbscan_clustering(['EV_sum', 'XLs_sum', 'Mem_score_sum', 'EM3D_EM'], min_cluster_size=10000)
AT.summarize_XLs_info()
exit()


