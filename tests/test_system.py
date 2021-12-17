import pyscal.core as pc
import os
import pyscal.crystal_structures as pcs
import numpy as np
from ase.build import bulk

def test_system_init():
	atoms, box = pcs.make_crystal(structure='bcc', 
                              lattice_constant=3.127, repetitions=(10,10,10),)
	sys = pc.System()
	sys.box = box
	sys.atoms = atoms

	assert len(sys.atoms["positions"]) == 10*10*10*2
	assert sys.triclinic == 0
	assert np.abs(sys.boxdims[0] - box[0][0]) < 1E-5


def test_system_triclinic():
	struct = bulk('Cu').repeat(10)
	sys = pc.System()
	sys.box = np.array(struct.cell)
	atoms = {}
	atoms["positions"] = struct.positions
	sys.atoms = atoms
	assert sys.triclinic == 1
	tb = (sys.rot== np.array(struct.cell).T)
	assert np.prod(tb) == 1
	tb = (sys.rotinv==np.linalg.inv(np.array(struct.cell).T))
	assert np.prod(tb) == 1

def test_nop():
	atoms, box = pcs.make_crystal(structure='bcc', 
                              lattice_constant=3.127, repetitions=(2,2,2),)
	sys = pc.System()
	sys.box = box
	sys.atoms = atoms

	assert sys.natoms == 16
	assert len(sys.atoms['positions']) == 2000
