import IMP
import IMP.pmi.output
import glob
import time
import os.path
import itertools
import sys


def extract_replica_info(run_folder):
"""
    function analyze the exhaustiveness of the replica exchange sampling
    % low exchange probability btx replicas means exchange was not sufficient
    If all runs have such case then rerun with smaller temperature difference between max_temp and min_temp 
    or increase number of replicas
    
    run_folder: path of the folder containing the output of your structural sampling
"""
# 1. extract the info form the stat files
    out_files = run_folder + "/output/stat_replica.*.out"
#    print(out_files)
    rex_out_files=glob.glob(out_files)
    temp_key="ReplicaExchange_CurrentTemp"
    maxtf_key="ReplicaExchange_MaxTempFrequency"
    mintf_key="ReplicaExchange_MinTempFrequency"
    ssr_key="ReplicaExchange_SwapSuccessRatio"
    score_key="score"

    # 2. calculate avg temperature of each replica
    score_temp_dict = {}
    avtemps_replicas=[]
    avg_swap_ratio_success = 0.0
    for f in rex_out_files:
        o=IMP.pmi.output.ProcessOutput(f)
        d=o.get_fields([temp_key,maxtf_key,mintf_key,
                    ssr_key,score_key])
        temps=[float(f) for f in d[temp_key]]
        scores=[float(f) for f in d[score_key]]
        suceess_ratio = d[ssr_key][-1]
        avg_swap_ratio_success += float(suceess_ratio)
        avtemp=sum(temps)/len(temps)
        avtemps_replicas.append(avtemp)
        print(f, avtemp)
        for n,t in enumerate(temps):
            s=scores[n]
            if t not in score_temp_dict:
                score_temp_dict[t]=[s]
            else:
                score_temp_dict[t].append(s)

    # 3. calculate average swapping success
    print("avg % swap success: ", avg_swap_ratio_success/len(rex_out_files))

    # 4. test that the average temperature per replica are similar
    # and check how probable the exchange between replica i with all other replica
    # if temperature difference is high, exchange probability is low
    counter = 0.0
    low_exchange = 0.0
    for c in itertools.combinations(avtemps_replicas,2):
        counter += 1
        if abs(c[0] - c[1]) > 0.1:
            low_exchange += 1
            
    # lower is better
    print(" % low exchange probability between replicas", (low_exchange/counter)*100)


# It is assumed that multiple independent conformational samplings were performed and all results of those runs are stored in folders named as run1, run2, etc ... 
# assuming all run folders are in your current directory (otherwise change the next line)

run_folder_path='./run*'
run_folders = glob.glob(run_folder_path)
for i in range(len(run_folders)):   
    print("extract info: ", run_folders[i])
    extract_replica_info(run_folders[i])
    print("----end-----")

sys.exit()
