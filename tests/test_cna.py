import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

#test conventional cna for fcc
def test_cna_cutoff():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.01)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    vec = sys.calculate_cna(cutoff=0.8536*4.1)
    assert vec[1] == 7*7*7*4

def test_cna_a1():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.01)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    vec = sys.calculate_cna(cutoff=None)
    assert vec[1] == 7*7*7*4

#now test adaptive
def test_cna_adaptive():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.calculate_cna(cutoff=None)
    atoms = sys.atoms
    assert atoms[0].structure == 1

    atoms, box = pcs.make_crystal(structure="hcp", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.calculate_cna(cutoff=None)
    atoms = sys.atoms
    assert atoms[0].structure == 2

    atoms, box = pcs.make_crystal(structure="bcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.calculate_cna(cutoff=None)
    atoms = sys.atoms
    assert atoms[0].structure == 3
