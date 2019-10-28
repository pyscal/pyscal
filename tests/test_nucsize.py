import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_system_nucsize():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    #test that atoms are set properly
    assert len(sys.atoms) == 16

    #now calculate nucsize
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    assert 16 == sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)

    #check the atom cluster props
    #assert np.round(sys.get_connection(atoms[0], atoms[1]), decimals=2) == 1.00

    #nothing to assert - just check if it works
    sys.calculate_frenkelnumbers()
    sys.find_clusters()
    assert sys.find_largest_cluster() == 16

def test_complex_system():
    sys = pc.System()
    sys.read_inputfile('examples/cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    assert 176 == sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)

def test_system_nucsize_fraction():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    #test that atoms are set properly
    assert len(sys.atoms) == 16

    #now calculate nucsize
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    assert 16 == sys.find_solids(bonds=0.8, threshold=0.5, avgthreshold=0.6, cluster=True)
