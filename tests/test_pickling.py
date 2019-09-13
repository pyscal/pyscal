import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import pyscal.pickle as pp


def test_pickle_atoms():

    sys = pc.System()
    sys.read_inputfile('tests/conf.dump')
    atoms = sys.get_atoms()

    #pickle and unpickel atom
    rn = np.random.randint(0, len(atoms)-1)
    patom = pp.pickle_atom(atoms[rn])
    assert patom.pos == atoms[rn].get_x()
    uatom = pp.unpickle_atom(patom)
    #assert uatom.get_x() == atoms[0].get_x()

    #pickle array of atoms
    pp.write_atoms('tests/pk.npy', atoms)
    ratoms = pp.read_atoms('tests/pk.npy')

    assert len(ratoms) == len(atoms)
    assert ratoms[rn].get_x() == atoms[rn].get_x()

    if os.path.exists('tests/pk.npy'):
        os.remove('tests/pk.npy')


def test_pickle_system():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [1, 1, 1])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    sys.find_neighbors(method = 'voronoi')

    #test write and read system
    sys.to_file('tests/sy.npy')

    #now read the pickled system
    psys = pc.System()
    psys.from_file('tests/sy.npy')

    #now get atoms and a random number of atom
    satoms = sys.get_atoms()
    patoms = psys.get_atoms()

    rn = np.random.randint(0, len(satoms)-1)
    assert satoms[rn].get_neighbors() == patoms[rn].get_neighbors()

    if os.path.exists('tests/sy.npy'):
        os.remove('tests/sy.npy')

    #now finally test pickling of series of atoms
    #syss = [sys, sys]
    #pp.write_systems('tests/sys.npy', syss)
    #psyss = pp.read_systems('tests/sys.npy')
    #assert len(psyss) == len(syss)

    #if os.path.exists('tests/sys.npy'):
    #    os.remove('tests/sys.npy')



