import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_triclinic():
    sys = pc.System()
    sys.read_inputfile('tests/conf.primitive.bcc.supercell.dump')
    sys.find_neighbors(method = 'cutoff', cutoff=1.2)
    atoms = sys.get_atoms()
    neighs = atoms[0].get_neighbors()
    assert len(neighs) == 14

    sys.find_neighbors(method = 'cutoff', cutoff=0.9)
    atoms = sys.get_atoms()
    neighs = atoms[0].get_neighbors()
    assert len(neighs) == 8
