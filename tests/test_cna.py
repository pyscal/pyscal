import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

#test conventional cna for fcc
def test_cna_cutoff():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.calculate_cna(cutoff=0.8536*4)
    atoms = sys.atoms
    assert atoms[0].cna[0] == [4, 2, 1, 12]

def test_cna_cutoff():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box
    sys.calculate_cnavector()
    atoms = sys.atoms
    assert atoms[0].cna[0] == [4, 2, 1, 12]

#now test adaptive
def test_cna_cutoff():
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
