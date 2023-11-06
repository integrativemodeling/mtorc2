#!/usr/bin/env python3

import numpy
import sys
import scipy.stats
import pandas as pd
import math

# Many of these functions are inspired from https://github.com/salilab/imp-sampcon

def get_contingency_table(cluster_info, total_num_models):
    num_clusters = len(cluster_info)
    full_ctable=numpy.zeros((num_clusters,2))

    cluster_members = []
    cluster_map ={}

    cluster_ranked = sorted(cluster_info, key=lambda k: len(cluster_info[k]), reverse=True)
    indx = 0
    for keys in cluster_info:
        cluster_map[str(indx)] =  keys
        cluster_members.append(cluster_info[keys])
        indx += 1

    for ic,cluster in enumerate(cluster_members):
        for member in cluster:
            if member.startswith('h1'):
                #print("run1", model_index)
                full_ctable[ic][0]+=1.0
            elif member.startswith('h2'):
                #print("run2", model_index)
                full_ctable[ic][1]+=1.0
            else:
                print('name of the members should start with h1 or h2')
                print('otherwise change this function')
                sys.exit()
    ## now normalize by number of models in each run
    numModelsRun1 = float(numpy.sum(full_ctable,axis=0)[0])
    numModelsRun2 = float(numpy.sum(full_ctable,axis=0)[1])
    print('total samples from run 1:',  int(numModelsRun1), 'total samples from run 2:', int(numModelsRun2))
    reduced_ctable=[]
    retained_clusters=[]

    for i in range(num_clusters):
        if full_ctable[i][0]<=10.0 or full_ctable[i][1]<=10.0:
        #if full_ctable[i][0]<=0.10*numModelsRun1 and full_ctable[i][1] <= 0.10*numModelsRun2:
            continue
        reduced_ctable.append([full_ctable[i][0],full_ctable[i][1]])
        retained_clusters.append(cluster_map[str(i)])
    clusters_with_num_models = {'cluster_name': [], 'percentage_of_all_members' : []}
    for i in range(0, len(retained_clusters)):
        clusters_with_num_models['percentage_of_all_members'].append(len(cluster_info[retained_clusters[i]])*100/total_num_models)
        clusters_with_num_models['cluster_name'].append(cluster_info[retained_clusters[i]][0])
    df = pd.DataFrame.from_dict(clusters_with_num_models).sort_values(by=['percentage_of_all_members'], ascending=False)
    print(df)
#    print(retained_clusters)
    return numpy.array(reduced_ctable),retained_clusters

def test_sampling_convergence(contingency_table,total_num_models):
    if len(contingency_table)==0:
        return 0.0,1.0

    ct = numpy.transpose(contingency_table)
    [chisquare,pvalue,dof,expected]=scipy.stats.chi2_contingency(ct)

#    print('difference in number of samples between expected and current_ctable \n', numpy.transpose(numpy.round(expected-ct)))
    if dof==0.0:
        cramersv=0.0
    else:
        cramersv=math.sqrt(chisquare/float(total_num_models))

    return(pvalue,cramersv)


def percent_ensemble_explained(ctable,total_num_models):
    if len(ctable)==0:
        return 0.0
    percent_clustered=float(numpy.sum(ctable,axis=0)[0]+numpy.sum(ctable,axis=0)[1])*100.0/float(total_num_models)
    return percent_clustered


def file_handler(filename):
    cluster_info ={}
    total_num_models = 0

    with open(filename, 'r') as read_file:
        for lines in read_file:
            if lines.startswith('cluster_number'):
                total_num_models += 1
                line_content = lines.strip().split(' ')
                cluster_id = 'cluster_' + line_content[1]
                first_member = line_content[3]
                cluster_info[cluster_id] = [first_member.split('/')[-1]]
            if 'cluster_mem' in lines:
                total_num_models += 1
                cluster_mem = lines.strip().split(' ')[1]
                cluster_info[cluster_id].append(cluster_mem.split('/')[-1])
    return total_num_models, cluster_info


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    print("provide a file name")


#filename = 'XX_cluster_info_details_cluster_radius_20.txt' # this file should be generated after running 03_structure_based_clustering.py
total_num_models, cluster_info = file_handler(filename)

#print(len(cluster_info))
c_table, retained_clusters = get_contingency_table(cluster_info, total_num_models)
pvalue, cramersv = test_sampling_convergence(c_table, total_num_models)
percent_clustered = percent_ensemble_explained(c_table,total_num_models)
# sampling_precision

if percent_clustered>80.0:
            if pvalue>0.05 or cramersv<0.10:
                print('sampling_precision is radius of the cluster')
            else:
                print('no success, use bigger cluster radius')

print('pvalue :', pvalue, 'Cramers V :', cramersv)

