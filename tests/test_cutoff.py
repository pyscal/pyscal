import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_cutoff():
    atoms, box = pcs.make_crystal(structure="fcc", lattice_constant=4.07, repetitions=(6,6,6))
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.find_neighbors(method="cutoff", cutoff=0)
    sys.set_atom_cutoff(factor=2)
    atoms = sys.atoms
    assert (atoms[0].cutoff - 5.755849198858498) < 1E-5
