#!/usr/bin/env python

from ast import literal_eval
import time

import IMP
import IMP.rmf
import RMF
import numpy as np
import glob
import sys
import os
import itertools
import shutil
import time


def align_models (model_1, model_2, trans_type, seq_ali, cluster_id):
    """
    A function to calculate the minimum rmsd between two models
    This models has C3 symmetry and translational symmetry along the C3 symmetry axis
    model_1: reference model
    model_2: query model
    seq_ali: a tuple of tuple. It represents way to aligning two models
    cluster_id
    """
    # model_1 and model_2 should be two rmf3 files

    fh_ref   = RMF.open_rmf_file_read_only(model_1)
    fh_query = RMF.open_rmf_file_read_only(model_2)
    mdl_ref  = IMP.Model()
    mdl_query  = IMP.Model()
    root_hier_ref = IMP.rmf.create_hierarchies(fh_ref, mdl_ref)[0]
    root_hier_query = IMP.rmf.create_hierarchies(fh_query, mdl_query)[0]

    # align the query model with respect to the reference model
    identifier = trans_type.split('_')
    if 'rot2' in identifier:
        chain_id_list = ['B', 'C', 'A']
    elif 'rot1' in identifier:
        chain_id_list = ['C', 'A', 'B']
    else:
        chain_id_list = ['A', 'B', 'C']
    # A look up Table
    table_ref = []
#    table_ref = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
    table_query = []
    for indx_1, indx_2 in seq_ali:
                    ref_chn_A = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][0]).get_selected_particles()
                    table_ref += ref_chn_A[8:] # not collecting first 8 particles, since they are part of a disordered region

                    ref_chn_B = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][1]).get_selected_particles()
                    table_ref += ref_chn_B[8:]

                    ref_chn_C = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][2]).get_selected_particles()
                    table_ref += ref_chn_C[8:]


                    query_chn_A = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][chain_id_list.index('A')]).get_selected_particles()
                    table_query += query_chn_A[8:]

                    query_chn_B = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][chain_id_list.index('B')]).get_selected_particles()
                    table_query += query_chn_B[8:]

                    query_chn_C = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][chain_id_list.index('C')]).get_selected_particles()
                    table_query += query_chn_C[8:]

    transformation = IMP.atom.get_transformation_aligning_first_to_second(table_query, table_ref)
    IMP.atom.transform(root_hier_query, transformation) # transform the hierarchy in place

    # create a new rmf3 file
    new_file = os.path.basename(model_2).split(".rmf3")[0] + "_aligned.rmf3"
    cwd = os.getcwd()
    path_new = cwd + "/ZZ_cluster" + cluster_id + "/"
    new_file_path = path_new + new_file
#    print(new_file_path)
    if not os.path.exists('ZZ_cluster' + cluster_id):
        os.makedirs('ZZ_cluster' + cluster_id)

    fh_new = RMF.create_rmf_file(new_file_path)
    IMP.rmf.add_hierarchy(fh_new, root_hier_query)

    IMP.rmf.save_frame(fh_new)
   

    del fh_ref
    del fh_query


# My models contain a Fibril with C3 symmetry and translational 
# There are 15 3-fold chains

chain_names_starts = ["A", "B", "C"] # 0 degrees rotation along C3 axis

chain_names = []

#assuming there are 15 chains and the chains are numbered as !!!!! CHECK YOUR .rmf3 FILES !!!!!

#[['A1', 'B1', 'C1'], ['A2', 'B2', 'C2'], ['A3', 'B3', 'C3'], ['A4', 'B4', 'C4'], ['A5', 'B5', 'C5'], ['A6', 'B6', 'C6'], ['A7', 'B7', 'C7'], ['A', 'B', 'C'], ['A8', 'B8', 'C8'], 
#['A9', 'B9', 'C9'], ['A10', 'B10', 'C10'], ['A11', 'B11', 'C11'], ['A12', 'B12', 'C12'], ['A13', 'B13', 'C13'], ['A14', 'B14', 'C14']]

for f in range(0, 15):
    chain_names.append([])
    if f < 7:
       for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e]+str(f+1))

    elif f >= 8:
        for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e]+str(f))

    else:
        for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e])


file_name = 'WWW_rmsf_info_10.txt'
transformation_dic = {}
head_list = ['h1_run28_38700.rmf3',  'h2_run12_3480.rmf3', 'h2_run50_12980.rmf3',  'h2_run52_8480.rmf3', 'h2_run27_44440.rmf3']
start = time.time()
with open(file_name, 'r') as read_file:
    for lines in read_file:
        line = lines.split('[')
        mem_info =  line[-2].split(',')
        rmsd = mem_info[0]
        cluster_head = mem_info[2].replace('"', '').split('/')[-1].split("'")[0]  
        if cluster_head in head_list:
            trans_type = line[0].split('(')[-1].strip("'").split("'")[0]
            cluster_mem = mem_info[1].replace('"', '').split('/')[-1].split("'")[0]
            transformation = literal_eval(line[-1].split(']')[0].strip('\''))
            if cluster_head in transformation_dic:
                transformation_dic[cluster_head].append([cluster_mem, trans_type, transformation, rmsd])
            else:
                transformation_dic[cluster_head] = [[cluster_mem, trans_type, transformation, rmsd]]
cluster_id = 0
for f in head_list:
    model_1 = f
    cluster_id += 1 
    for k in transformation_dic[f]:
        model_2 = k[0]
        trans_type = k[1]
        transformation = k[2]
        start = time.time()
        align_models (model_1, model_2, trans_type, transformation, str(cluster_id))
        end = time.time()
        print("time required", end-start)
#print(transformation_dic['h1_run28_38700.rmf3'][0])
