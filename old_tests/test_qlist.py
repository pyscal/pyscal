import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_q_list():
    #this might take a while, it will find all qs
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    sys.find_neighbors(method = 'voronoi')

    sys.calculate_q([2, 4], averaged=True)
    q = sys.get_qvals([2, 4], averaged=True)
    #x1 = q[0]
    #x2 = q[1]
    #xx1 = np.round(np.mean(np.array(x1)), decimals=2)
    #xx2 = np.round(np.mean(np.array(x2)), decimals=2)
    #assert xx1 == 0.00
    #assert xx2 == 0.22
