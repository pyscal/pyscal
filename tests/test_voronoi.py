import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_voro_props():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)

    sys.find_neighbors(method = 'voronoi')
    atoms = sys.get_atoms()
    atom = atoms[0]
    v = atom.get_vorovector()
    assert v == [0,6,0,8]