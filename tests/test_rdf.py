import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_rdf():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims
    rdf, r = sys.calculate_rdf()

    args = np.argsort(rdf)
    assert(np.round(r[args][-1], decimals=2) == 0.87)

    #test for an fcc type
    atoms, boxdims = pcs.make_crystal('fcc', repetitions = [4, 4, 4])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims
    rdf, r = sys.calculate_rdf()

    args = np.argsort(rdf)
    assert(np.round(r[args][-1], decimals=2) == 0.69)
