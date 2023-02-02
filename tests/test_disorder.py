import pytest
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_ordered_disorder():
    sys = pc.System('tests/files/conf.fcc.dump')
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_q(6)
    sys.calculate_disorder(averaged=True)
    assert np.mean(sys.disorder) < 0.50
    assert np.mean(sys.avg_disorder) < 0.50
    assert np.mean(sys.atom.steinhardt.disorder.norm) < 0.50
    assert np.mean(sys.atom.steinhardt.disorder.average) < 0.50
    
def test_disordered_disorder():
    sys = pc.System('tests/files/conf.lqd.dump')
    sys.find_neighbors(method='cutoff', cutoff=0)
    sys.calculate_q(6)
    sys.calculate_disorder(averaged=True)
    assert np.mean(sys.disorder) > 1.00
    assert np.mean(sys.avg_disorder) > 1.00
