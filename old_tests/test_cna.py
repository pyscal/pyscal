import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
from ase.build import bulk

#test conventional cna for fcc
def test_cna_cutoff():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.01)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    res = sys.calculate_cna(lattice_constant=4.00)
    assert res["fcc"] == 7*7*7*4

def test_cna_a1():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.01)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    res = sys.calculate_cna(lattice_constant=None)
    assert res["fcc"] == 7*7*7*4

#now test adaptive
def test_cna_adaptive():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.calculate_cna(lattice_constant=None)
    atoms = sys.atoms
    assert atoms[0].structure == 1

    atoms, box = pcs.make_crystal(structure="hcp", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.calculate_cna(lattice_constant=None)
    atoms = sys.atoms
    assert atoms[0].structure == 2

    atoms, box = pcs.make_crystal(structure="bcc", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.calculate_cna(lattice_constant=None)
    atoms = sys.atoms
    assert atoms[0].structure == 3

def test_cna_diamond():
    atoms, box = pcs.make_crystal(structure="diamond", repetitions=(7,7,7), lattice_constant=4.00, noise=0.1)
    sys = pc.System()
    sys.box = box
    sys.atoms = atoms
    sys.identify_diamond()
    atoms = sys.atoms
    assert atoms[0].structure == 5


def test_ase_bulks():
    
    al_fcc = bulk("Al")
    fe_bcc = bulk("Fe")
    ti_hcp = bulk("Ti")

    sys = pc.System()
    sys.read_inputfile(al_fcc, format="ase")
    cna = sys.calculate_cna()
    assert cna["fcc"] == 1

    sys = pc.System()
    sys.read_inputfile(fe_bcc, format="ase")
    cna = sys.calculate_cna()
    assert cna["bcc"] == 1

    sys = pc.System()
    sys.read_inputfile(ti_hcp, format="ase")
    cna = sys.calculate_cna()
    assert cna["hcp"] == 2