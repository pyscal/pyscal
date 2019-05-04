import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs

def test_system_nucsize():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)

    #test that atoms are set properly
    assert len(sys.get_allatoms()) == 16

    #now calculate nucsize
    sys.set_nucsize_parameters(0.867, 6, 0.5, 0.5)
    #sys.calculate_nucsize()
    assert sys.calculate_nucsize() == 16

    #check the atom cluster props
    atoms = sys.get_allatoms()
    atom = atoms[0]
    cluster = atom.get_cluster()
    assert atom.get_cluster() == [1,0,1,1]
    assert sys.get_largestcluster() == 1
    #assert np.round(sys.get_connection(atoms[0], atoms[1]), decimals=2) == 1.00

    #nothing to assert - just check if it works
    sys.calculate_frenkelnumbers()
    sys.find_clusters()
    assert sys.find_largest_cluster() == 16