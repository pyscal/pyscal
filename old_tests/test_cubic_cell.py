import pytest
import os
import numpy as np
import pyscal.core as pc
from ase.build import bulk

def test_cubic():
    #this might take a while, it will find all qs
    atoms = bulk("Fe")
    sys = pc.System()
    sys.read_inputfile(atoms, format="ase")
    box, atoms = sys.extract_cubic_box()

    assert np.sum(np.array(atoms[0].pos)-np.array([0.0, 0.0, 0.0])) < 1E-10
    assert np.sum(np.array(atoms[1].pos)-np.array([1.435, 1.435, 1.435])) < 1E-10
