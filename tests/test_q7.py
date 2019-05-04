import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs

def test_q_7():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    sys.get_neighbors(method = 'voronoi')

    sys.calculate_q(7, averaged=True)
    q = sys.get_qvals(7, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00  