import pytest
import numpy as np
import pyscal.core as pc
from pyscal.crystal_structures import Structure

def test_chiparamsbcc():
    sys = Structure().lattice.bcc(repetitions = [3, 3, 3], lattice_constant=4)
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip2 = [3, 0, 0, 0, 36, 12, 0, 36, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0


def test_chiparamsfcc():
    sys = Structure().lattice.fcc(repetitions = [5, 5, 5], lattice_constant=4)
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip2 =  [6, 0, 0, 0, 24, 12, 0, 24, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0

def test_chiparamsdia():
    sys = Structure().lattice.diamond(repetitions = [3, 3, 3], lattice_constant=4)
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip2 =  [0, 0, 0, 0, 6, 0, 0, 0, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0
