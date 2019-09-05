import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_q_3():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    #sys.get_neighbors(method = 'voronoi')
    sys.get_neighbors(method = 'cutoff', cutoff=0.9)
    
    sys.calculate_q(3, averaged=True)
    q = sys.get_qvals(3, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00    