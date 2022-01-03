import pytest
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_chiparamsbcc():
    atoms, box = pcs.make_crystal('bcc', repetitions=[3,3,3], lattice_constant=4)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip = sys.chiparams[2]
    chip2 = [3, 0, 0, 0, 36, 12, 0, 36, 0]
    assert np.sum(np.array(chip)-np.array(chip2)) == 0


def test_chiparamsfcc():
    atoms, box = pcs.make_crystal('fcc', repetitions=[5,5,5], lattice_constant=4)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip = sys.chiparams[2]
    chip2 =  [6, 0, 0, 0, 24, 12, 0, 24, 0]
    assert np.sum(np.array(chip)-np.array(chip2)) == 0

def test_chiparamsdia():
    atoms, box = pcs.make_crystal('diamond', repetitions=[3,3,3], lattice_constant=4)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_chiparams()
    chip = sys.chiparams[2]
    chip2 =  [0, 0, 0, 0, 6, 0, 0, 0, 0]
    assert np.sum(np.array(chip)-np.array(chip2)) == 0
