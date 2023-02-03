import pytest
import os
import numpy as np
import pyscal.core as pc
from pyscal.crystal_structures import Structure

def test_q_4():
    sys = Structure().lattice.bcc(repetitions = [4, 4, 4])

    #sys.get_neighbors(method = 'voronoi')
    sys.find_neighbors(method = 'cutoff', cutoff=0.9)

    q = sys.calculate_q([4, 6])
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51 , "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63 , "Calculated q4 value is wrong!"

    q = sys.calculate_q([4, 6], averaged=True)
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51 , "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63 , "Calculated q4 value is wrong!"

