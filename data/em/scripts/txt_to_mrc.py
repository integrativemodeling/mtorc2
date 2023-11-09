import IMP
import IMP.atom
import IMP.isd.gmm_tools
import sys,os

input_txt = sys.argv[1]
ref_mrc = sys.argv[2]
out_prefix = input_txt.split('.txt')[0]

print(input_txt, ref_mrc)
# Read gmms from txt file

mdl = IMP.Model()
ps = []
gmms = IMP.isd.gmm_tools.decorate_gmm_from_text(input_txt,
                                                ps,
                                                mdl)

p = ps[0]
print(IMP.core.Gaussian(p).get_gaussian())
print(IMP.atom.Mass(p).get_mass())
shape = IMP.core.Gaussian(p).get_gaussian()
covar=[c for row in IMP.algebra.get_covariance(shape) for c in row]
mean=list(shape.get_center())

print(gmms, type(ps[0]), covar, mean, shape.get_variances())



# Get bounding box
ref_map = IMP.em.read_map(ref_mrc,IMP.em.MRCReaderWriter())
bb = IMP.em.get_bounding_box(ref_map)

# Write mrc
IMP.isd.gmm_tools.write_gmm_to_map(ps, out_prefix+'.mrc',1.668)

