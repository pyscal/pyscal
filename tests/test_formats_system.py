import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_lammps_dump():
    sys = pc.System()
    sys.read_inputfile('tests/conf.dump')
    atoms = sys.atoms
    assert len(atoms) == 500

    #check box
    assert sys.box == [[-7.66608, 11.1901],[-7.66915, 11.1931],[-7.74357, 11.2676]]

    #check few atoms
    filtered_atoms = [ atom for atom in atoms if atom.id == 204]
    assert filtered_atoms[0].pos == [-0.10301, -6.35752, -6.44787]

    #now check the same for zipped file
    sys = pc.System()
    sys.read_inputfile('tests/conf.dump.gz')
    atoms = sys.atoms
    assert len(atoms) == 500

    #check box
    assert sys.box == [[-7.66608, 11.1901],[-7.66915, 11.1931],[-7.74357, 11.2676]]

    #check few atoms
    filtered_atoms = [ atom for atom in atoms if atom.id == 204]
    assert filtered_atoms[0].pos == [-0.10301, -6.35752, -6.44787]

def test_neighbors_number():
    atoms, box = pcs.make_crystal(structure="fcc", repetitions=(3,3,3), lattice_constant=4.00)
    sys = pc.System()
    sys.atoms = atoms
    sys.box = box

    sys.find_neighbors(method="number", nmax=8)
    atoms = sys.atoms
    nns = np.mean([atom.coordination for atom in atoms])
    assert nns == 8

    sys.find_neighbors(method="number", nmax=12)
    atoms = sys.atoms
    nns = np.mean([atom.coordination for atom in atoms])
    assert nns == 12


def test_scaled():
    sys = pc.System()
    sys.read_inputfile('tests/conf.bcc.scaled.dump')
    atoms = sys.atoms
    assert len(atoms) == 2
    assert atoms[1].pos[0] == 1.0

def test_poscar():
    sys = pc.System()
    sys.read_inputfile('tests/POSCAR', format='poscar')
    atoms = sys.atoms
    assert len(atoms) == 42

    #now assert atoms of different types
    type1 = len([atom for atom in atoms if atom.type == 1])
    type2 = len([atom for atom in atoms if atom.type == 2])
    type3 = len([atom for atom in atoms if atom.type == 3])
    assert type1 == 38
    assert type2 == 2
    assert type3 == 2

    #now test the coorfinates of the atom
    assert sys.box == [[0.0, 18.768662916], [0.0, 8.9728430088], [0.0, 2.83746]]

    #now check coordinates of first atom
    assert atoms[0].pos == [0.020389021322710754, 8.8981229427339, 0.0005978263145216028]
