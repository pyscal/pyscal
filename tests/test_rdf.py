import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_rdf_bcc():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [10, 10, 10])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms
    rdf, r = sys.calculate_rdf(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.86 < 1E-5)

def test_rdf_fcc():
    atoms, boxdims = pcs.make_crystal('fcc', repetitions = [10, 10, 10])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms
    rdf, r = sys.calculate_rdf(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.70 < 1E-5)
