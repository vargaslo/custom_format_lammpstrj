import numpy as np
import itertools

def subtraj_ind(s):

    NF, Nmolec = np.shape(s)

    num_regions = 2  # micro region is 1, meso region is 2)
    region = np.ones_like(s)  # default to region 1

    for i in range(num_regions-1):
        radius = 17  # mesochannel has 17 A radius
        region += (s <= radius)  # increase region number if closer to corner


#    myregions = {'micro': 1, 'meso': 2}
    myregions = {1: 'micro', 2: 'meso'}

    # get indices that bound the subtrajectories
    indices = {}
    indices['micro'] = []
    indices['meso'] = []
    first = 0
    last = 0
    for p in range(Nmolec):
        for v in itertools.groupby(region[:,p]):
            key, group = v
            N = len(list(group))
            last = first + N

            indices[myregions[int(key)]].append((first, last))
            first = last

    # convert bounding indices to actual indices
    hashes = {}
    hashes['micro'] = []
    hashes['meso'] = []
    for k,v in indices.items():
        for i in v:
            all_indices = np.arange(i[0], i[1]).tolist()
            hashes[k].append(all_indices)

    return hashes



def assemble_list_of_trjs():

    list_of_trjs = []
