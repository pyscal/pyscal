import pyscal.core as pc
import os
from pyscal.crystal_structures import Structure
import numpy as np
from ase.build import bulk

def test_system_init():
    sys = Structure().lattice.bcc(repetitions = [10,10,10], lattice_constant=3.127)
    sys.find_neighbors(method="cutoff", cutoff=3.6)
    a1 = np.array(sys.atoms.neighbors.distance)
    a2 = np.array([2.708061437633939,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 2.708061437633937,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339356])
    assert np.sum(a1-a2) < 1E-5
    assert np.sum(sys.atoms.neighbors.distance[0]-a2) < 1E-5

    sys.find_neighbors(method="cutoff", cutoff=3.6)
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939,
 3.126999999999999,
 3.126999999999999,
 3.126999999999999,
 3.127,
 3.127,
 3.127])
    assert np.sum(a1-a2) < 1E-5

    sys.find_neighbors(method="cutoff", cutoff='sann')
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939,
 3.126999999999999,
 3.126999999999999,
 3.126999999999999,
 3.127,
 3.127,
 3.127,
 4.422245809540667])
    assert np.sum(a1-a2) < 1E-5

    sys.find_neighbors(method="number", nmax=8)
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939])
    assert np.sum(a1-a2) < 1E-5