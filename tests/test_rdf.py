import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
from pyscal.crystal_structures import Structure


def test_rdf_bcc():
    sys = Structure().lattice.bcc(repetitions = [10, 10, 10])
    rdf, r = sys.calculate_rdf(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.86 < 1E-5)

def test_rdf_fcc():
    sys = Structure().lattice.fcc(repetitions = [10, 10, 10])
    rdf, r = sys.calculate_rdf(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.70 < 1E-5)
