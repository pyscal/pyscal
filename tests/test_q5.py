import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_q_5():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    sys.find_neighbors(method = 'voronoi')

    sys.calculate_q(5, averaged=True)
    q = sys.get_qvals(5, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00
