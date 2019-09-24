import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import pyscal.pickle as pp



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
    assert satoms[rn].neighbors == patoms[rn].neighbors

    if os.path.exists('tests/sy.npy'):
        os.remove('tests/sy.npy')

    #now finally test pickling of series of atoms
    #syss = [sys, sys]
    #pp.write_systems('tests/sys.npy', syss)
    #psyss = pp.read_systems('tests/sys.npy')
    #assert len(psyss) == len(syss)

    #if os.path.exists('tests/sys.npy'):
    #    os.remove('tests/sys.npy')
