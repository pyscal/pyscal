import pytest
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_ordered_disorder():
    sys = pc.System()
    sys.read_inputfile('examples/conf.fcc.dump')
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_q(6)
    sys.calculate_disorder(averaged=True)
    atoms = sys.atoms
    disorder = [atom.disorder for atom in atoms]
    assert np.mean(disorder) < 0.50

    disorder = [atom.avg_disorder for atom in atoms]
    assert np.mean(disorder) < 0.50

def test_disordered_disorder():
    sys = pc.System()
    sys.read_inputfile('examples/conf.lqd')
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_q(6)
    sys.calculate_disorder(averaged=True)
    atoms = sys.atoms
    disorder = [atom.disorder for atom in atoms]
    assert np.mean(disorder) > 1.00

    disorder = [atom.avg_disorder for atom in atoms]
    assert np.mean(disorder) > 1.00
