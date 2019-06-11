import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs

def test_q_2():
    #this might take a while, it will find all qs
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    #sys.read_inputfile("tests/bcc.dat")
    sys.get_neighbors(method = 'voronoi')
    atoms = sys.get_allatoms()
    vols = []
    avgvols = []
    for atom in atoms:
        vols.append(atom.get_volume())
        avgvols.append(atom.get_volume(averaged=True))
    assert np.mean(np.array(vols)) == 0.5
    assert np.mean(np.array(avgvols)) == 0.5