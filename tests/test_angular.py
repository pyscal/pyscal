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

    sys.find_neighbors(method = 'cutoff', cutoff=0)
    sys.calculate_angularcriteria()

    assert np.round(np.mean(np.array(sys.atoms.angular_parameters.diamond_angle)), decimals=2) == 0.00
