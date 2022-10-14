import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_angular():
    atoms, boxdims = pcs.make_crystal('diamond', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    sys.find_neighbors(method = 'voronoi')
    sys.calculate_angularcriteria()

    q = [atom.angular for atom in sys.atoms]
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00
