import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_q_4():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims
    #sys.get_neighbors(method = 'voronoi')
    sys.find_neighbors(method = 'cutoff', cutoff=0.9)
    sys.calculate_q([4, 6], averaged=True)
    atoms = sys.atoms
    #assert np.round(atoms[0].q[2], decimals=2) == 0.51
    assert np.round(atoms[0].get_q(4, averaged=True), decimals=2) == 0.51
    assert np.round(atoms[0].get_q([4, 6])[0], decimals=2) == 0.51
    assert np.round(atoms[0].get_q([4, 6], averaged=True)[1], decimals=2) == 0.63
    assert np.round(atoms[0].get_q([4, 6]), decimals=2)[0] == 0.51
    assert np.round(atoms[0].get_q([4, 6], averaged=True)[1], decimals=2) == 0.63

    #now change the q4 values
    atoms[0].set_q(4, .23)
    atoms[0].set_q(4, .23, averaged=True)
    assert np.round(atoms[0].get_q(4), decimals=2) == 0.23
    assert np.round(atoms[0].get_q(4, averaged=True), decimals=2) == 0.23

    atoms[0].set_q([4, 6], [.23, .46])
    atoms[0].set_q([4, 6], [.23, .46], averaged=True)
    assert np.round(atoms[0].get_q(4), decimals=2) == 0.23
    assert np.round(atoms[0].get_q(4, averaged=True), decimals=2) == 0.23
    assert np.round(atoms[0].get_q(6), decimals=2) == 0.46
    assert np.round(atoms[0].get_q(6, averaged=True), decimals=2) == 0.46
