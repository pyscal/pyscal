import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_cs_ges():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4], lattice_constant=4.00)
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    sys.calculate_centrosymmetry(nmax=8)
    atoms = sys.atoms
    q = [atom.centrosymmetry for atom in atoms]
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00

    sys.calculate_centrosymmetry(nmax=12)
    atoms = sys.atoms
    q = [atom.centrosymmetry for atom in atoms]
    assert np.round(np.mean(np.array(q)), decimals=2) > 0.00

def test_cs_gvm():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [4, 4, 4], lattice_constant=4.00)
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    sys.calculate_centrosymmetry(nmax=8, algorithm="gvm")
    atoms = sys.atoms
    q = [atom.centrosymmetry for atom in atoms]
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00

    sys.calculate_centrosymmetry(nmax=12, algorithm="gvm")
    atoms = sys.atoms
    q = [atom.centrosymmetry for atom in atoms]
    assert np.round(np.mean(np.array(q)), decimals=2) > 0.00

def test_csm_ovito():
    sys = pc.System()
    sys.read_inputfile("tests/bcc.csm.dump", customkeys=["Centrosymmetry"])    
    sys.calculate_centrosymmetry(nmax=8)
    atoms = sys.atoms
    ind = np.random.randint(0, len(atoms))
    assert np.abs(atoms[ind].centrosymmetry - float(atoms[ind].custom["Centrosymmetry"])) < 1E-5
