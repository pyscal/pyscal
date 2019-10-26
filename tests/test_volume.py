import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_q_2():
    #this might take a while, it will find all qs
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    #sys.read_inputfile("tests/bcc.dat")
    sys.find_neighbors(method = 'voronoi')
    atoms = sys.atoms
    vols = []
    avgvols = []
    for atom in atoms:
        vols.append(atom.volume)
        avgvols.append(atom.avg_volume)
    assert np.mean(np.array(vols)) == 0.5
    assert np.mean(np.array(avgvols)) == 0.5
