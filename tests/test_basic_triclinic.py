import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_triclinic():
    sys = pc.System()
    sys.read_inputfile('tests/conf.primitive.bcc.supercell.dump')
    sys.find_neighbors(method = 'cutoff', cutoff=1.2)
    atoms = sys.atoms
    neighs = atoms[0].neighbors
    assert len(neighs) == 14

    sys.find_neighbors(method = 'cutoff', cutoff=0.9)
    atoms = sys.atoms
    neighs = atoms[0].neighbors
    assert len(neighs) == 8

def test_triclinic_frames():

    os.system("cat tests/conf.primitive.bcc.supercell.dump tests/conf.primitive.bcc.supercell.dump > tests/bcc.prim.dat")
    sys3 = pc.System()
    sys3.read_inputfile("tests/bcc.prim.dat")
    sys3.find_neighbors(method = 'cutoff', cutoff=0.9)
    atoms = sys3.atoms
    neighs = atoms[0].neighbors
    assert len(neighs) == 8
