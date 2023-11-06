import os
import shutil
import numpy as np
from IMP import ArgumentParser
import scipy as sp
from scipy import spatial

import IMP.sampcon
from IMP.sampcon import scores_convergence, clustering_rmsd
from IMP.sampcon import rmsd_calculation, precision_rmsd

import IMP
import sys
import glob
import RMF
import sys



def get_rmfs_coordinates(str_file, chain_names, subunit_name, num=0):

#    conform = []
#    num = 0
    masses = []
    radii = []
    ps_names = []

    models_name = []

#    print(str_file)
    models_name.append(str_file)

    m = IMP.Model()
    inf = RMF.open_rmf_file_read_only(str_file)
    h = IMP.rmf.create_hierarchies(inf, m)[0]
    IMP.rmf.load_frame(inf, 0)

    pts = []

    s0 =[]
    if subunit_name == 'Abeta':
        for indx_1 in range(len(chain_names)):
            for indx_2 in range(len(chain_names[0])):
                s0 += IMP.atom.Selection(h, resolution=1,molecule=subunit_name, chain_id = chain_names[indx_1][indx_2]).get_selected_particles()[8:]
    else:
        s0 = IMP.atom.Selection(h, resolution=1, molecule=subunit_name).get_selected_particles()
#        s0 = IMP.atom.Selection(h, resolution=1).get_selected_particles()
#    print(s0)
    for leaf in s0:

        p=IMP.core.XYZR(leaf)
        pts.append([p.get_coordinates()[i] for i in range(3)])

        if num == 0:
            masses.append(IMP.atom.Mass(leaf).get_mass())
            radii.append(p.get_radius())
        mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(leaf))
        # traverse up the Hierarchy to get copy number of the molecule that the bead belongs to
        copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(leaf))

        if IMP.atom.Fragment.get_is_setup(leaf): #TODO not tested on non-fragment systems
            residues_in_bead = IMP.atom.Fragment(leaf).get_residue_indexes()

            ps_names.append(mol_name+"_"+str(min(residues_in_bead))+"_"+str(max(residues_in_bead))+"_"+str(copy_number))

        else:
            residue_in_bead = str(IMP.atom.Residue(leaf).get_index())

            ps_names.append(mol_name+"_"+residue_in_bead+"_"+residue_in_bead+"_"+str(copy_number))

    conform = pts
    pts = []
    num = num + 1

    return ps_names, masses, radii, np.array(conform), models_name



# My models contain a Fibril with C3 symmetry and translational 
# There are 15 3-fold chains

chain_names_starts = ["A", "B", "C"] # 0 degrees rotation along C3 axis
chain_names = [[], [], []]

#assuming there are 15 chains and the chains are numbered as !!!!! CHECK YOUR .rmf3 FILES !!!!!

#[['A1', 'B1', 'C1'], ['A2', 'B2', 'C2'], ['A3', 'B3', 'C3'], ['A4', 'B4', 'C4'], ['A5', 'B5', 'C5'], ['A6', 'B6', 'C6'], ['A7', 'B7', 'C7'], ['A', 'B', 'C'], ['A8', 'B8', 'C8'], 
#['A9', 'B9', 'C9'], ['A10', 'B10', 'C10'], ['A11', 'B11', 'C11'], ['A12', 'B12', 'C12'], ['A13', 'B13', 'C13'], ['A14', 'B14', 'C14']]
for e in range(len(chain_names_starts)):
    for f in range(0, 15):
        if f < 7:
                chain_names[e].append(chain_names_starts[e]+str(f+1))

        elif f >= 8:
                chain_names[e].append(chain_names_starts[e]+str(f))

        else:
                chain_names[e].append(chain_names_starts[e])

# Obtain the subunits for which we need to calculate densities
density = 'density.txt'
density_custom_ranges = precision_rmsd.parse_custom_ranges(density)

# somr parameters

density_threshold = 30
voxel = 3
#align = True
align = None
symmetry_groups = False
# List of cluster members
i = 0
cluster_head = 'h1_run28_38700.rmf3' # change this to some other filename # this is cluster center of Cluster I 
cluster_member_list = []
ist_rmf = sorted(glob.glob("*.rmf3"))
for f in ist_rmf:
    if f != cluster_head:
        cluster_member_list.append(f)

# Output cluster precisions
fpc = open("Cluster_Precision.txt", 'w+')


# The cluster centroid is the first conformation.
# We use this as to align and compute RMSD/precision
ps_names, masses, radii, conform_0, models_name = get_rmfs_coordinates(cluster_head, chain_names, subunit_name='TREM2', num=0)
#print(ps_names)
# create a directory for the cluster
if not os.path.exists("./cluster.%s" %i):
    os.mkdir("./cluster.%s" %i)
    os.mkdir("./cluster.%s/Sample_A/" % i)
    os.mkdir("./cluster.%s/Sample_B/" % i)
else:
    shutil.rmtree("./cluster.%s" %i)
    os.mkdir("./cluster.%s" %i)
    os.mkdir("./cluster.%s/Sample_A/" % i)
    os.mkdir("./cluster.%s/Sample_B/" % i)

gmd1 = precision_rmsd.GetModelDensity(
            custom_ranges=density_custom_ranges,
            resolution=density_threshold, voxel=voxel,
            bead_names=ps_names)
gmd2 = precision_rmsd.GetModelDensity(
            custom_ranges=density_custom_ranges,
            resolution=density_threshold, voxel=voxel,
            bead_names=ps_names)
gmdt = precision_rmsd.GetModelDensity(
        custom_ranges=density_custom_ranges,
        resolution=density_threshold, voxel=voxel,
        bead_names=ps_names)

# Also output the identities of cluster members
both_file = open('cluster.'+str(i)+'.all.txt', 'w')
sampleA_file = open('cluster.'+str(i)+'.sample_A.txt', 'w')
sampleB_file = open('cluster.'+str(i)+'.sample_B.txt', 'w')

# Create a model with just the cluster_member particles
model = IMP.Model()
ps = [] # particle list to be updated by each RMF frame
for pi in range(len(conform_0)):
    p = IMP.Particle(model, "%s" % str(pi))
    IMP.core.XYZ.setup_particle(p, (0,0,0))
    IMP.core.XYZR.setup_particle(p, float(radii[pi]))
    IMP.atom.Mass.setup_particle(p, float(masses[pi]))
    ps.append(p)

# Obtain cluster precision by obtaining average RMSD of each model
# to the cluster center
cluster_precision = 0.0


# symmetry group (defining as None):

# transformation from internal pyRMSD orientation
trans = None
# for each model in the cluster (modify)
for mem in cluster_member_list:
    # mem is a rmf file
    # conform_0 coordinates of the cluster head
    ps_names_, masses_, radii_, conform_mem, models_name_ = get_rmfs_coordinates(mem, chain_names, subunit_name='TREM2', num=1)
#    model_index = all_models[mem]
    # get superposition of each model to cluster center and the
    # RMSD between the two
    if symmetry_groups:
        rmsd, superposed_ps, trans = \
            precision_rmsd.get_particles_from_superposed_amb(
                conform_mem, conform_0, align, ps, trans,
                symm_groups)
    else:
         rmsd, superposed_ps, trans = \
            precision_rmsd.get_particles_from_superposed(
                conform_mem, conform_0, align, ps, trans)
    model.update() # why not?

    cluster_precision += rmsd
    gmdt.add_subunits_density(superposed_ps) #superposed_ps # total density map
    print(mem, file=both_file)

    if mem.startswith('h1'):
        # density map for sample A
        gmd1.add_subunits_density(superposed_ps)
        print(mem, file=sampleA_file)
    else:
        # density map for sample B
        gmd2.add_subunits_density(superposed_ps)
        print(mem, file=sampleB_file)

cluster_precision /= float(len(cluster_member_list))

print("Cluster precision (average distance to cluster centroid) "
            "of cluster with cluster head ", cluster_head, " is %.3f" % cluster_precision, "A",
              file=fpc)
both_file.close()
sampleA_file.close()
sampleB_file.close()

# Output density files for the cluster
density = gmdt.write_mrc(path="./cluster.%s" % i, file_prefix = "LPD")
gmd1.write_mrc(path="./cluster.%s/Sample_A/" % i, file_prefix = "LPD")
gmd2.write_mrc(path="./cluster.%s/Sample_B/" % i, file_prefix = "LPD")

fpc.close()

