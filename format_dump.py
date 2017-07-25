#
#  File expects THREE command line arguments
#
#  the first command line argument should be a string
#  that points to the dump file
#
#  the second command line argument is the time step (fs per timestep)
#
#  the third command line argument is the desired number of frames to
#  retrieve from the dump file

import sys
import os
import read_dump
import numpy as np
import coord_transform
import itertools
import subtrajectories

#------------------------------------------------------------------------------

# dumpfile should be first command line argument
dumpfile = os.path.abspath(sys.argv[1])

# timestep is user-defined!
timestep_fs = float(sys.argv[2])

# desired number of frames (can be 'All' or an integer)
try:
    desired_NF = int(sys.argv[3])
except:
    desired_NF = sys.argv[3]

# crystal parameters for NU-1000
crystal = {}
crystal['a'] = 39.97
crystal['b'] = 40.00
crystal['c'] = 16.58*2
crystal['alpha'] = 90
crystal['beta'] = 90
crystal['gamma'] = 120

#------------------------------------------------------------------------------

# check if dumpfile exists and is a file
if not os.path.isfile(dumpfile):
    raise FileNotFoundError
else:
    dirpath = os.path.dirname(dumpfile)

# read in dumpfile and save to dict
out = read_dump.read_dumpfile(dumpfile, timestep_fs=timestep_fs, NF=desired_NF)

# assign variable names to dict items
t_ps = out['t_ps']
delta_t = t_ps[1] - t_ps[0]  # ps per frame
xyz_unwrapped = out['xyz']
Nmolec = out['mlc']
NF = out['numframes']

# use crystal parameters to find transformation matrices
Mfwd = coord_transform.xyz2fracM(crystal)
Mrev = coord_transform.frac2xyzM(crystal)

# use transformation matrices to find wrapped coordinates
xyz_wrapped = coord_transform.wrapcoords(xyz_unwrapped, Mfwd, Mrev)
s = coord_transform.xy2s(xyz_wrapped, Mrev)

# export (composite) trajectories to file
outfile = '{}.composite'.format(dumpfile)
with open(outfile, 'w') as fout:
    header = '# trj xw yw zw s Zu  ---  ps per frame: {}\n'.format(delta_t)
    fout.write(header)

    for trj in range(Nmolec):
        for frame in range(NF):
            xw, yw, zw = xyz_wrapped[frame, trj, :]
            Xu, Yu, Zu = xyz_unwrapped[frame, trj, :]
            s_dist = s[frame, trj]
            fout.write('{:d} '.format(trj))
            fout.write('{:.4f} '.format(xw))
            fout.write('{:.4f} '.format(yw))
            fout.write('{:.4f} '.format(zw))
            fout.write('{:.4f} '.format(s_dist))
            fout.write('{:.4f} '.format(Zu))
            fout.write('\n')


if 1:  # custom script for NU-1000 to split micro/meso regions

    # use s parameter to create subtrajectories
    subtrajs = subtrajectories.subtraj_ind(s)

    # export (micro- and meso- ) subtrajectories to file
    for k, v in subtrajs.items():
        outfile = '{}.{}'.format(dumpfile, k)
        with open(outfile, 'w') as fout:
            header = '# trj xw yw zw s Zu  ---  ps per frame: {}\n'.format(delta_t)
            fout.write(header)

            trj_i = 0
            for subtraj in v:
                c, r = np.unravel_index(subtraj, (NF, Nmolec), order='F')
                xw, yw, zw = np.transpose(xyz_wrapped[c, r, :])
                Xu, Yu, Zu = np.transpose(xyz_unwrapped[c, r, :])
                s_dist = s[c, r]
                for i in range(len(subtraj)):
                    fout.write('{:d} '.format(trj_i))
                    fout.write('{:.4f} '.format(xw[i]))
                    fout.write('{:.4f} '.format(yw[i]))
                    fout.write('{:.4f} '.format(zw[i]))
                    fout.write('{:.4f} '.format(s_dist[i]))
                    fout.write('{:.4f} '.format(Zu[i]))
                    fout.write('\n')
                trj_i += 1
