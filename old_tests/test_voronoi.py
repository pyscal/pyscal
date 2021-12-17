import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_voro_props():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [10, 10, 10])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms


    sys.find_neighbors(method = 'voronoi')
    sys.calculate_vorovector()
    atoms = sys.atoms
    atom = atoms[0]
    #v = atom.vorovector
    #assert v == [0,6,0,8]
