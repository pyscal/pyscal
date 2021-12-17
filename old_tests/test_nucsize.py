import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_complex_system():
    sys = pc.System()
    sys.read_inputfile('examples/cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    assert 176 == sys.find_solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)

def test_cluster():
    sys = pc.System()
    sys.read_inputfile('examples/cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    sys.find_solids(cluster=False)
    val = sys.cluster_atoms("solid", largest = True)
    assert 176 == val

def test_cluster_cutoff():
    sys = pc.System()
    sys.read_inputfile('examples/cluster.dump')
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    sys.find_solids(cluster=False)
    val = sys.cluster_atoms("solid", largest = True, cutoff=3.63)
    assert 176 == val

def test_system_nucsize_fraction():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2], lattice_constant=3.20)
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    #test that atoms are set properly
    assert len(sys.atoms) == 16

    #now calculate nucsize
    sys.find_neighbors(method='cutoff', cutoff=3.63)
    assert 16 == sys.find_solids(bonds=0.8, threshold=0.5, avgthreshold=0.6, cluster=True)