import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs


def test_file_system():

    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [1, 1, 1])
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms
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
    sys3.box = boxdims
    sys3.atoms = atoms
    sys3.to_file('tests/tjkf.dat', customkeys=['velocity'])

    #now read it again
    sys4 = pc.System()
    sys4.read_inputfile('tests/tjkf.dat', customkeys=['velocity'])
    #now get atoms and check them
    atoms = sys4.atoms
    assert int(atoms[0].custom['velocity']) == 12
    assert int(atoms[1].custom['velocity']) == 24
    sys4.to_file('tests/tjkf.dat.gz', customkeys=['velocity'],
        compressed=True)

    #now read it again
    sys5 = pc.System()
    sys5.read_inputfile('tests/tjkf.dat.gz', customkeys=['velocity'])
    #now get atoms and check them
    atoms = sys5.atoms
    assert int(atoms[0].custom['velocity']) == 12
    assert int(atoms[1].custom['velocity']) == 24

    sys6 = pc.System()
    sys6.read_inputfile('tests/tjkf.dat')
    sys6.to_file('tests/poscar_1', format="poscar", species=['Mo'])
    sys6.to_file('tests/data_1', format="lammps-data", species=['Mo'])