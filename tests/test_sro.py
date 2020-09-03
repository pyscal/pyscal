import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_sro():
    atoms, box = pcs.make_crystal('l12', lattice_constant=4.00, repetitions=[2,2,2])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.find_neighbors(method='cutoff', cutoff=4.5)
    sro = sys.calculate_sro(reference_type=1, average=True)
    assert np.round(sro[0], decimals=2) == -0.33
    assert sro[1] == 1.0

    atoms = sys.atoms
    sro = atoms[4].sro
    assert np.round(sro[0], decimals=2) == -0.33
    assert sro[1] == 1.0

    sys.find_neighbors(method='cutoff', cutoff=4.5)
    sro = sys.calculate_sro(reference_type=1, average=True, shells=1)
    assert np.round(sro[0], decimals=2) == -0.07
