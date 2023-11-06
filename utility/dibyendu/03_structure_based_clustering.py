from __future__ import print_function
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
import argparse

parser = argparse.ArgumentParser(description='Cluster models based on a given cluster radius')
parser.add_argument('cluster_radius', metavar='cluster_radius', type=float,
                    help='a float representing the radius of a structural cluster provide in A unit')

args = parser.parse_args()

def get_coordinates(imp_particle_list): # we might not need this
    """
    extract coordinates from a particle list
    """
    coord_xyz = [IMP.core.XYZ(l) for l in imp_particle_list]
    return coord_xyz


def extract_particles_and_do_rmsd(align_list, root_hier_ref, root_hier_query, table_ref, table_query, chain_order, first_model=None):
    """
    Extract particles from certain chains of molecule = Abeta and calculate the rmsd between ref can query model
    align_list = which chain of query model to align with which chain of the reference model
    root_hier_ref = root hierarchy of the reference model
    root_hier_query = root hierarchy of the query model
    table_ref = A look up list for reference model 
    table_query = A look up list for query model
    chain_order = defines rotation symmetry ex. [A, B, C] no rotation, [B, C, A] 120 degrees rotation, [C, A, B] 240 degrees rotation
    first_model = name of the query (rmf3 file)
    """

    
    rmsd_ex = 0         # initialize a variable to store rmsd
    ref_par_list = []   # two lists to contain particles these lists will be used to
    query_par_list = [] # get coordinate of the particles using 'get_coordinates' function


    #align_list has a form: [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10), (11, 11), (12, 12), (13, 13), (14, 14)]
    # each tuple consists of two indices (refrence, query)
    

    for indx_1, indx_2 in align_list:

#        ref_chn_A = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = ref_chain_id_list[indx_1][0]).get_selected_particles()
        ref_par_list += table_ref[indx_1][0]
        ref_par_list += table_ref[indx_1][1]
        ref_par_list += table_ref[indx_1][2]

        
#        query_chn_A = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = query_chain_id_list[indx_2][0]).get_selected_particles()
        query_par_list += table_query[indx_2][chain_order.index('A')]
        query_par_list += table_query[indx_2][chain_order.index('B')]
        query_par_list += table_query[indx_2][chain_order.index('C')]

    # using get_coordinates function to get coordinates of a list of particles (Note: we might not need this step)
#    ref_all = get_coordinates(ref_par_list)
#    qu_all = get_coordinates(query_par_list)

    ref_TREM2 = IMP.atom.Selection(root_hier_ref, molecule='TREM2').get_selected_particles()

    query_TREM2 = IMP.atom.Selection(root_hier_query, molecule='TREM2').get_selected_particles()  # delete later

#    rmsd_ini = IMP.atom.get_rmsd(ref_TREM2, query_TREM2)

    # align first
    # 1. get the transformatiion vector

    transformation = IMP.atom.get_transformation_aligning_first_to_second(query_par_list, ref_par_list)

    IMP.atom.transform(root_hier_query, transformation) # transform the hierarchy in place

    # calculate rmsd

    rmsd_ex = IMP.atom.get_rmsd(ref_TREM2, query_TREM2)

    reversed_transform = transformation.get_inverse()

    IMP.atom.transform(root_hier_query, reversed_transform) # transform the hierarchy in place


    # align first molecule to the second. Here first molecule is query
    # get the transformation vector
    if first_model != None:
        transformation = IMP.atom.get_transformation_aligning_first_to_second(query_par_list, ref_par_list)
        IMP.atom.transform(root_hier_query, transformation) # transform the hierarchy in place

        # create a new rmf3 file
        new_file = os.path.basename(first_model).split(".rmf3")[0] + "_new.rmf3"
        cwd = os.getcwd()
        path_new = cwd + "/new/"
        new_file_path = path_new + new_file
        print(new_file_path)
        if not os.path.exists("new"):
            os.makedirs("new")

        fh_new = RMF.create_rmf_file(new_file_path)
        IMP.rmf.add_hierarchy(fh_new, root_hier_query)

        IMP.rmf.save_frame(fh_new)
    return rmsd_ex 

def shifted_rmsd_calculation (model_1, model_2, seq_ali, rmsd_threshold):
    """
    A function to calculate the minimum rmsd between two models
    This models has C3 symmetry and translational symmetry along the C3 symmetry axis
    model_1: reference model
    model_2: query model
    seq_ali: a dictionary of list. It represents all the difference ways aligning two models
    rmsd_threshold: a threshold to stop further alignment testing
    """
    # model_1 and model_2 should be two rmf3 files

    fh_ref   = RMF.open_rmf_file_read_only(model_1)
    fh_query = RMF.open_rmf_file_read_only(model_2)
    mdl_ref  = IMP.Model()
    mdl_query  = IMP.Model()
    root_hier_ref = IMP.rmf.create_hierarchies(fh_ref, mdl_ref)[0]
    root_hier_query = IMP.rmf.create_hierarchies(fh_query, mdl_query)[0]

    # A look up Table
    table_ref = []
#    table_ref = [[], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
    table_query = []
    for keys in seq_ali:
        if keys.startswith('forward_trans'):
                for indx_1, indx_2 in seq_ali['forward_trans'][0]:
                    table_ref.append([])
                    table_query.append([])
                    ref_chn_A = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][0]).get_selected_particles()
                    table_ref[indx_1].append(ref_chn_A[8:]) # not collecting first 8 particles, since they are part of a disordered region

                    ref_chn_B = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][1]).get_selected_particles()
                    table_ref[indx_1].append(ref_chn_B[8:])

                    ref_chn_C = IMP.atom.Selection(root_hier_ref, molecule="Abeta", chain_id = chain_names[indx_1][2]).get_selected_particles()
                    table_ref[indx_1].append(ref_chn_C[8:])


                    query_chn_A = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][0]).get_selected_particles()
                    table_query[indx_2].append(query_chn_A[8:])

                    query_chn_B = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][1]).get_selected_particles()
                    table_query[indx_2].append(query_chn_B[8:])

                    query_chn_C = IMP.atom.Selection(root_hier_query, molecule="Abeta", chain_id = chain_names[indx_2][2]).get_selected_particles()
                    table_query[indx_2].append(query_chn_C[8:])

#    print(table_query[0])


    rmsd_dic = {} # a dictionary to store rmsd for each alignment
    for keys in seq_ali:
        count = 0
        if keys.endswith("trans"):
            chain_order = ['A', 'B', 'C']
            for align_list in seq_ali[keys]:
                dic_keys = keys + str(count)
                count = count + 1
#                print(dic_keys)
                rmsd_ex = extract_particles_and_do_rmsd(align_list, root_hier_ref, root_hier_query, table_ref, table_query, chain_order, first_model=None)
                rmsd_dic[dic_keys] = [rmsd_ex, model_2, model_1, align_list]
            temp_rmsd_tuple = sorted(rmsd_dic.items(), key=lambda x:x[1][0])[0]
            temp_rmsd = round(temp_rmsd_tuple[1][0])
            if temp_rmsd < rmsd_threshold:
    #            final_rmsd = temp_rmsd_tuple
                break
        elif keys.endswith("rot1"):
            chain_order = ['C', 'A', 'B']
            count_rot1 = 0
            for align_list in seq_ali[keys]:
                rmsd_ex = 0.0
                dic_keys_rot1 = keys + "_" + str(count_rot1)
                count_rot1 = count_rot1 + 1
                rmsd_ex = extract_particles_and_do_rmsd(align_list, root_hier_ref, root_hier_query, table_ref, table_query, chain_order, first_model=None)
                rmsd_dic[dic_keys_rot1] = [rmsd_ex, model_2, model_1, align_list]
            
            temp_rmsd_tuple = sorted(rmsd_dic.items(), key=lambda x:x[1][0])[0]
            temp_rmsd = round(temp_rmsd_tuple[1][0])
            if temp_rmsd < rmsd_threshold:
    #            final_rmsd = temp_rmsd_tuple
                break
        elif keys.endswith("rot2"):
            chain_order = ['B', 'C', 'A']
            count_rot2 = 0
            for align_list in seq_ali[keys]:
                rmsd_ex = 0.0
                dic_keys_rot2 = keys + "_" + str(count_rot2)
                count_rot2 = count_rot2 + 1
                rmsd_ex = extract_particles_and_do_rmsd(align_list, root_hier_ref, root_hier_query, table_ref, table_query, chain_order, first_model=None)
                rmsd_dic[dic_keys_rot2] = [rmsd_ex, model_2, model_1, align_list]
    final_rmsd = sorted(rmsd_dic.items(), key=lambda x:x[1][0])[0]

    # align the query model with respect to the reference model
    """
    identifier = final_rmsd[0].split('_')[1]
    if identifier.endswith("rot2"):
        chain_id_list = ['B', 'C', 'A']
    elif identifier.endswith("rot1"):
        chain_id_list = ['C', 'A', 'B']
    else:
        chain_id_list = ['A', 'B', 'C']
    final_align_list = final_rmsd[1][3]
    fh_query4 = RMF.open_rmf_file_read_only(model_2)
    mdl_query4  = IMP.Model()
    root_hier_query4 = IMP.rmf.create_hierarchies(fh_query4, mdl_query4)[0]
    extract_particles_and_do_rmsd(final_align_list, root_hier_ref, root_hier_query, table_ref, table_query, chain_id_list, first_model=model_2)
    
    del mdl_query4
    del root_hier_query4
    """

    del fh_ref
    del fh_query

    return final_rmsd


def check_old_clusters(cluster_dic, current_id, mdl_2, write_file):
    """
    To Check if current model 2 could be part of any old clusters from the cluster_dic
    cluster_dic: dictionary of the exsisting clusters
    current_id: The cluster we do not want to check
    write_file: write transformation info
    """
    found_cluster = False
    if len(cluster_dic) !=0:
      d_new = sorted(cluster_dic, key=lambda k: len(cluster_dic[k]), reverse=True)
      for k in range(len(d_new)):
           keys = d_new[k]
           if keys != current_id:
               cluster_head = cluster_dic[keys][0]
               temp_rmsd_tuple = shifted_rmsd_calculation (cluster_head, mdl_2, seq_ali, rmsd_threshold)
#               print( "checking cluster_", keys)
               print('rmsd_info: ', temp_rmsd_tuple, file=write_file)
               temp_rmsd = round(temp_rmsd_tuple[1][0])
               if temp_rmsd < rmsd_threshold_cluster:
                  print("found cluster_",keys)
                  cluster_name = "cluster_"+ str(keys)
                  found_cluster = True
                  cluster_dic[keys].append(mdl_2)
                  break
    return found_cluster



#####THIS SECTION IS UNIQUE FOR MY PROJECT  ###########
#####ALWAYS CHECK THE CODE BETWEEN NEXT $$$ ###########

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# My models contain a Fibril with C3 symmetry and translational 
# There are 15 3-fold chains

chain_names_starts = ["A", "B", "C"] # 0 degrees rotation along C3 axis
rotated_1 = ["B", "C", "A"] # 120 degrees rotation 
rotated_2 = ["C", "A", "B"] # 240 degrees rotation

chain_names = []
rot1_chain_names = []
rot2_chain_names = []

#assuming there are 15 chains and the chains are numbered as !!!!! CHECK YOUR .rmf3 FILES !!!!!

#[['A1', 'B1', 'C1'], ['A2', 'B2', 'C2'], ['A3', 'B3', 'C3'], ['A4', 'B4', 'C4'], ['A5', 'B5', 'C5'], ['A6', 'B6', 'C6'], ['A7', 'B7', 'C7'], ['A', 'B', 'C'], ['A8', 'B8', 'C8'], 
#['A9', 'B9', 'C9'], ['A10', 'B10', 'C10'], ['A11', 'B11', 'C11'], ['A12', 'B12', 'C12'], ['A13', 'B13', 'C13'], ['A14', 'B14', 'C14']]

for f in range(0, 15):   
    chain_names.append([])
    rot1_chain_names.append([])
    rot2_chain_names.append([])
    if f < 7:
       for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e]+str(f+1))
           rot1_chain_names[f].append(rotated_1[e]+str(f+1))
           rot2_chain_names[f].append(rotated_2[e]+str(f+1))

    elif f >= 8:
        for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e]+str(f))
           rot1_chain_names[f].append(rotated_1[e]+str(f))
           rot2_chain_names[f].append(rotated_2[e]+str(f))

    else:
        for e in range(len(chain_names_starts)):
           chain_names[f].append(chain_names_starts[e])
           rot1_chain_names[f].append(rotated_1[e])
           rot2_chain_names[f].append(rotated_2[e])


#print(chain_names)
#print(rot1_chain_names)
#print(rot2_chain_names)

# dictionary to store indices

seq_ali = {'forward_trans': [], 'backward_trans': [], 'forward_rot1': [], 'backward_rot1': [], 'forward_rot2': [], 'backward_rot2': []}

for i in range(len(chain_names)):
    k = 0
    m = len(chain_names)-1
    demo_list1 = []
    for j  in range(i, len(chain_names)):
        demo_list1.append((j, k))
        k = k+1
        m = m-1
    seq_ali['forward_trans'].append(demo_list1)   
    seq_ali['forward_rot1'].append(demo_list1)   
    seq_ali['forward_rot2'].append(demo_list1)   
for i in range(len(chain_names), 0, -1):
    m = 0
    demo_list2 = []
    for j  in range(i, len(chain_names)):
        demo_list2.append((m, j))
        m = m+1
    if len(demo_list2) !=0:
       seq_ali['backward_trans'].append(demo_list2)
       seq_ali['backward_rot1'].append(demo_list2)
       seq_ali['backward_rot2'].append(demo_list2)
#print(seq_ali)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


print(sys.argv)
models = glob.glob("./*.rmf3")
curnt_dir = os.getcwd()

rmsd_threshold_cluster = 10.0

# selecting a rmsd_threshold to decide if we can stop searching for more alignments, the value should be less than equal to rmsd_threshold_cluster

rmsd_threshold = 10.0


if len(sys.argv) > 1:
    rmsd_threshold_cluster = args.cluster_radius
else:
    print("No cluster threshold provided, using 10 Angs")

i = 0

# for future transformation store the best alignment with rmsd values
transformation_write_file = open('ZZ_cluster_transformation_details_' + str(int(rmsd_threshold_cluster))+ '.txt', 'w')

#print('total rmf3 file: ', len(models))
#print('cluster radius :', rmsd_threshold_cluster)

#print('rmsd_info: ', '1. alignment_type', '2. rmsd', '3. query_model', '4. reference_model', '5. alignment_details')
transformation_write_file.writelines('rmsd_info: \t 1. alignment_type \t 2. rmsd \t 3. query_model \t 4. reference_model \t 5. alignment_details\n')
cluster_dic = {}
total_models = len(models)
while i <len(models)-1:
##    print(i)
    model_1 = models[i] # initialize model 1 (reference model)
    currnt_id = i
    cluster_dic[i] = [model_1] # defining cluster dictionary to keep infromation about files
    cluster_dir = "cluster_"+ str(i) 
    

    noise_level = 0 # to monitor how many times a bad reference model will be used
    
    model_2 = models[i+1:]
    
#    print('models left to check ', len(model_2))
    for mdl_2 in model_2:
        i = i + 1 # checking on more model
#        print("current model number: {} of total models {} ".format(str(i), str(total_models)))
##        start = time.time()
        temp_rmsd_tuple = shifted_rmsd_calculation (model_1, mdl_2, seq_ali, rmsd_threshold)  # calculating RMSD
        print('rmsd_info: ', temp_rmsd_tuple, file=transformation_write_file)
        temp_rmsd = round(temp_rmsd_tuple[1][0])
 
        if temp_rmsd > rmsd_threshold_cluster: # cheking, if we want to put in the current cluster
#            print("can not be part of cluster_", currnt_id)
        
            if check_old_clusters(cluster_dic, currnt_id, mdl_2, transformation_write_file): # checking if we can put in some old cluster
            
                if noise_level < 5: # If we can use the current cluster, no point of using the cluster head as reference 
#                    print("noise level still low") 
                    noise_level = noise_level + 1 # we will try current cluster head as reference one less time than 5
                else:
#                    print("breaking due to high noise level")
#                    print("Dont use the current model 1 need a new reference")
                    best_cluster_id = sorted(cluster_dic, key=lambda k: len(cluster_dic[k]), reverse=True)[0]
#                    print("using again the cluster", best_cluster_id, " as main model 1 cluster")
                    currnt_id = best_cluster_id
                    model_1 = cluster_dic[best_cluster_id][0]
                    noise_level = 0
            else:
#                print("can not be part of any existing cluster")
                cluster_dir = "cluster_"+ str(i)
           
                break
        else:
                cluster_dir = "cluster_"+ str(currnt_id)
#                print("copying to", cluster_dir)
                cluster_dic[currnt_id].append(mdl_2)

transformation_write_file.close()

# write cluster information
w_file = open('XX_cluster_info_details_' + str(int(rmsd_threshold_cluster))+ '.txt', 'w')
w_file.write(' '.join(sys.argv) + ' <---- cluster radius ' +'\n')
for keys in cluster_dic:
    w_file.write('cluster_number: ' + str(keys) + ' cluster_head: ' + cluster_dic[keys][0]+'\n')
    for i in range(len(cluster_dic[keys])):
        w_file.write('\t' + '\t' + ' cluster_mem: ' + cluster_dic[keys][i]+'\n')

w_file.close()
##    print("cluster heads", cluster_dic[keys][0])
sys.exit()

