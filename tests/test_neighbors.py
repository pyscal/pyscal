import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_neighbors_system():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)

    #test that atoms are set properly
    assert len(sys.get_atoms()) == 432

    #then lets find neighbors
    #cutoff method - first shell only
    sys.find_neighbors(method = 'cutoff', cutoff=0.867)
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 8

    sys.reset_neighbors()

    #cutoff method - second shell
    sys.find_neighbors(method = 'cutoff', cutoff=1.1)
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 14

    #cutoff method - third shell
    sys.find_neighbors(method = 'cutoff', cutoff=1.5)
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 26

    sys.reset_neighbors()
    #voronoi method - first shell only
    sys.find_neighbors(method = 'voronoi')
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 14

    assert np.round(sys.get_distance(atoms[0], atoms[1]), decimals=2) == 0.87


def test_neighbors_system_filter():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [2, 2, 2])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)


    sys.find_neighbors(method = 'cutoff', cutoff=0.867)
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 8

    #now we take all the neighbors of atom0 and replace half of
    #them with a different type
    neighs = atoms[0].neighbors

    #replace the neighs
    atoms[neighs[0]].type = 2
    atoms[neighs[1]].type = 2

    #now set these atoms back
    #for atom in atoms:
    sys.set_atom(atoms[neighs[0]])
    sys.set_atom(atoms[neighs[1]])

    #recalculate neighbors with filter
    sys.find_neighbors(method = 'cutoff', cutoff=0.867, filter='type')
    #any atom should have 8 neighbors
    atoms = sys.get_atoms()
    assert atoms[0].coordination == 6
