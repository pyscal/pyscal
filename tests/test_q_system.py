import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_q_4():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    #sys.get_neighbors(method = 'voronoi')
    sys.get_neighbors(method = 'cutoff', cutoff=0.9)
    sys.calculate_q([4, 6], averaged=True)
    q = sys.get_qvals([4, 6], averaged=True)
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51 , "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63 , "Calculated q4 value is wrong!"   

    q = sys.get_qvals([4, 6])
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51 , "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63 , "Calculated q4 value is wrong!"   
