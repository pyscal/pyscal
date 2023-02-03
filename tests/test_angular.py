import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
from pyscal.crystal_structures import Structure

def test_angular():
    sys = Structure().lattice.diamond(repetitions = [4, 4, 4])
    sys.find_neighbors(method = 'cutoff', cutoff=0)
    sys.calculate_angularcriteria()

    assert np.round(np.mean(np.array(sys.atoms.angular_parameters.diamond_angle)), decimals=2) == 0.00
