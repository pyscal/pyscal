import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs

def test_neighbors_system():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)

    #then lets find neighbors
    #cutoff method - first shell only
    sys.get_neighbors(method = 'cutoff', cutoff=0)
    #any atom should have 8 neighbors
    atoms = sys.get_allatoms()
    assert atoms[0].get_coordination() == 14

    sys.get_neighbors(method = 'cutoff', cutoff=0, threshold=1)
    #any atom should have 8 neighbors
    atoms = sys.get_allatoms()
    assert atoms[0].get_coordination() == 14

    sys.calculate_q(8, averaged=True)
    q = sys.get_qvals(8, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.45


    sys.calculate_q(7, averaged=True)
    q = sys.get_qvals(7, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.05     