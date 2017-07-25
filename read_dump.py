#
#    Read LAMMPS dump file line by line
#    Check that molecule order is correct
#    Check that file is not corrupted
#
#    output dict contains:
#    out['timestep_per_frame'] is a user-defined time step used to calculate out['t_ps']
#    out['t_ps'] is a 1-d array with the time points
#    out['numframes'] is NF, the number of frames read in from the file
#    out['mlc'] is Nmol, the number of molecules
#    out['xyz'] is an array with dimesions given by (NF, Nmol, dim)
#


import os
import numbers
import collections
import numpy as np

def read_dumpfile(infile, timestep_fs=None, NF=None):

    # require user to input the timestep
    if not isinstance(timestep_fs, numbers.Real):
        raise TypeError('Timestep must be entered as argument to read_file function')

    t = []
    xyz = []
    tmp_xyz = []
    frame = -1

    with open(infile, 'rb') as fin:
        for line in fin:
            line_chunks = line.split()
            if line_chunks[0] != '#':

                if len(line_chunks)==4:
                    ind, x, y, z = map(float, line_chunks)

                    # check if molecules are ordered as expected
                    if ind==expected_N.popleft():
                        tmp_xyz.append([x, y, z])
                        # after reading all N molecules, append to array and reset tmp

                        if ind==N:
                            xyz.append(tmp_xyz)
                            t.append(timestep)
                            tmp_xyz = []

                            # stop reading if desired number of frames have been read
                            if isinstance(NF, numbers.Real) and len(t)>=NF:
                                break
                    else:
                        print ('\nWARNING: Molecule ordering seems off. File may be corrupted')
                        break

                if len(line_chunks)==2:

                    # use these to get delta timesteps
                    if frame==0:
                        timestep_0 = timestep
                    if frame==1:
                        delta_timestep = timestep - timestep_0

                    # guess what the next timestep should be
                    if frame>0:
                        expected_timestep = timestep + delta_timestep

                    # read new timestep and number of molec
                    timestep, N = [int(i) for i in line_chunks]

                    # check if timestep is as expected. otherwise break
                    if frame>0:
                         if timestep!=expected_timestep:
                             print ('WARNING: Partial file read. Frames are inconsistent.\n')
                             break

                    # expected number of trjs to read
                    expected_N = collections.deque([i+1 for i in range(N)])

                    frame+=1

        # Convert timestep to time
        fs_per_timestep = timestep_fs
        ps_per_fs = 0.001
        t_ps = (np.array(t) - t[0]) * fs_per_timestep * ps_per_fs

        # Convert to ndarray and get number of frames read
        xyz = np.array(xyz)
        numframes, N, dims = np.shape(xyz)

        print ('{} frames read for {} molecules from {}\n'.format(numframes, N, os.path.abspath(infile)))

    out = {}
    out['xyz'] = xyz
    out['numframes'] = numframes
    out['mlc'] = N
    out['t_ps'] = t_ps
    out['timesteps_per_frame'] = delta_timestep

    return out
