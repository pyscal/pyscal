import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_voro_props():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims


    sys.find_neighbors(method = 'voronoi')
    sys.calculate_vorovector()
    atoms = sys.atoms
    atom = atoms[0]
    v = atom.vorovector
    assert v == [0,6,0,8]
