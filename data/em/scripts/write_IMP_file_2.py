import numpy as np
import pandas as pd
import sys

input = sys.argv[1]

cols = ['g', 'W', 'M0', 'M1', 'M2', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz']
GMMs = pd.DataFrame(columns = cols)

gaussian = np.nan
g = -1
for line in open(input,'r'):
    if line[0:6] == 'HETATM':
        gaussian = int(line[6:12])
    elif line[0:6] == 'REMARK' and line[7:12] == 'GAUSS':
        vals = line.split()
        g = int(vals[2])
        if 'W' in vals:
            W = vals[4]
        elif 'M' in vals:
            M0 = vals[4]
            M1 = vals[5]
            M2 = vals[6]
        elif 'CovM' in vals and 'xx' in vals:
            xx = float(vals[5])
            xy = float(vals[7])
            xz = float(vals[9])
        elif 'CovM' in vals and 'yy' in vals:
            yy = float(vals[5])
            yz = float(vals[7])
            zz = float(vals[9])
    if g != -1 and gaussian != g:
        # GMMs = GMMs.append(pd.Series([g, W, M0, M1, M2, xx, xy, xz, yy, yz, zz], index=cols), ignore_index= True)
        GMMs.loc[len(GMMs)] = g, W, M0, M1, M2, xx, xy, xz, yy, yz, zz

        del W, M0, xx, xy, xz, yy, yz, zz
        g = -1


print(GMMs.head())
output = input.split('.gmm')[0]+'.txt'
out = open(output,'w')
out.write('#|num|weight|mean|covariance matrix| \n')
for i, row in GMMs.iterrows():
    num = row['g']
    w = row['W']
    m0 = row['M0']
    m1 = row['M1']
    m2 = row['M2']
    xx = row['xx']
    xy = row['xy']
    xz = row['xz']
    yy = row['yy']
    yz = row['yz']
    zz = row['zz']

    out.write(f'|{num}|{w}|{m0} {m1} {m2}|{xx} {xy} {xz} {xy} {yy} {yz} {xz} {yz} {zz}|\n')

out.close()

