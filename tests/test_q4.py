import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_q_4():
    atoms, boxdims = pcs.make_crystal('fcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    #sys.get_neighbors(method = 'voronoi')
    sys.find_neighbors(method = 'cutoff', cutoff=0.9)
    sys.calculate_q(4)
    q = sys.get_qvals(4)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.19 , "Calculated q4 value is wrong!"

    sys.calculate_q(4, only_averaged=True)
    q = sys.get_qvals(4, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.19 , "Calculated q4 value is wrong!"

def test_q_4_voro():
    atoms, boxdims = pcs.make_crystal('fcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    sys.find_neighbors(method = 'voronoi')
    sys.calculate_q(4)
    q = sys.get_qvals(4)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.19 , "Calculated q4 value is wrong!"

    sys.calculate_q(4, only_averaged=True)
    q = sys.get_qvals(4, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.19 , "Calculated q4 value is wrong!"
