import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import pyscal.pickle_object as pp



def test_pickle_system():

    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [1, 1, 1])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims
    sys.find_neighbors(method = 'voronoi')

    #test write and read system
    sys.to_pickle('tests/sy.npy')

    #now read the pickled system
    psys = pc.System()
    psys.from_pickle('tests/sy.npy')

    #now get atoms and a random number of atom
    satoms = sys.atoms
    patoms = psys.atoms

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


def test_file_system():

    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [1, 1, 1])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims
    sys.find_neighbors(method = 'voronoi')

    sys.to_file('tests/tjkf.dat')

    #now try to read in the file
    sys2 = pc.System()
    sys2.read_inputfile('tests/tjkf.dat')
    assert len(sys2.atoms) == 2

    #now add some custom values
    atoms[0].custom = {"velocity":12}
    atoms[1].custom = {"velocity":24}

    #now try to read in the file
    sys3 = pc.System()
    sys3.atoms = atoms
    sys3.box = boxdims
    sys3.to_file('tests/tjkf.dat', custom=['velocity'])

    #now read it again
    sys4 = pc.System()
    sys4.read_inputfile('tests/tjkf.dat', customkeys=['velocity'])
    #now get atoms and check them
    atoms = sys4.atoms
    assert int(atoms[0].custom['velocity']) == 12
    assert int(atoms[1].custom['velocity']) == 24
