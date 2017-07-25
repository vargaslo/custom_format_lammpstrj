import numpy as np

def cosd(deg):
    return np.cos(deg/180. * np.pi)

def sind(deg):
    return np.sin(deg/180. * np.pi)

def xyz2fracM(crystal_prms):
    # use crystal parameters to compute the transition matrix Mfwd
    # that will convert Cartesian coords to fractional coords
    a = crystal_prms['a']
    b = crystal_prms['b']
    c = crystal_prms['c']
    alpha = crystal_prms['alpha']
    beta = crystal_prms['beta']
    gamma = crystal_prms['gamma']
    v = np.sqrt(1-cosd(alpha)**2-cosd(beta)**2-cosd(gamma)**2 + 2*cosd(alpha)* cosd(beta)*cosd(gamma))
    r1 = [1./a, -cosd(gamma)/a/sind(gamma), (cosd(alpha)*cosd(gamma)-cosd(beta)) / a/v/sind(gamma)]
    r2 = [0, 1./b/sind(gamma), (cosd(beta)*cosd(gamma)-cosd(alpha)) / b/v/sind(gamma)]
    r3 = [0, 0, sind(gamma)/c/v]
    Mfwd = np.array([r1, r2, r3])
    return Mfwd

def frac2xyzM(crystal_prms):
    # use crystal parameters to compute the transition matrix Mrev
    # that will convert fractional coords to Cartesian coords
    a = crystal_prms['a']
    b = crystal_prms['b']
    c = crystal_prms['c']
    alpha = crystal_prms['alpha']
    beta = crystal_prms['beta']
    gamma = crystal_prms['gamma']
    v = np.sqrt(1-cosd(alpha)**2-cosd(beta)**2-cosd(gamma)**2 + 2*cosd(alpha)* cosd(beta)*cosd(gamma))
    r1 = [a, b*cosd(gamma), c*cosd(beta)]
    r2 = [0, b*sind(gamma), c*(cosd(alpha)-cosd(beta)*cosd(gamma))/sind(gamma)]
    r3 = [0, 0, c*v/sind(gamma)]
    Mrev = np.array([r1, r2, r3])
    return Mrev

def wrapcoords(alldata, Mfwd, Mrev):
    N, trj, dim = np.shape(alldata)
    alldataT= np.transpose(alldata)
    xyz_wrap_T = np.empty((dim, trj, N))
    for i in range(trj):
        frac_T = np.dot(Mfwd, alldataT[:,i,:])

        # wrap unit cell coordinates
        frac_wrap_T = frac_T % 1

        # convert fractional to cartesian coords
        xyz_wrap_T[:, i, :] = np.dot(Mrev, frac_wrap_T)

    xyz_wrap = np.transpose(xyz_wrap_T)

    return xyz_wrap

def xy2s(wrapped, Mrev):

    xyz = wrapped

    # get Cartesian coordinates of unit cell corners
    P1 = np.dot(Mrev, np.array([0.0, 0.0, 0]))   #   P4------P3
    P2 = np.dot(Mrev, np.array([1.0, 0.0, 0]))   #    \        \
    P3 = np.dot(Mrev, np.array([1.0, 1.0, 0]))   #     \   uc   \
    P4 = np.dot(Mrev, np.array([0.0, 1.0, 0]))   #      \        \
                                                      #       P1-------P2

    # find the distance closest to a corner
    s1 = np.sqrt(np.sum((xyz[:,:,0:2]-P1[np.newaxis, np.newaxis, 0:2])**2, axis=2))
    s2 = np.sqrt(np.sum((xyz[:,:,0:2]-P2[np.newaxis, np.newaxis, 0:2])**2, axis=2))
    s3 = np.sqrt(np.sum((xyz[:,:,0:2]-P3[np.newaxis, np.newaxis, 0:2])**2, axis=2))
    s4 = np.sqrt(np.sum((xyz[:,:,0:2]-P4[np.newaxis, np.newaxis, 0:2])**2, axis=2))
    s0 = np.minimum(s1, s2)
    s0 = np.minimum(s0, s3)
    s0 = np.minimum(s0, s4)

    return s0
