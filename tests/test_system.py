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

	for a in sys.iter_atoms():
		assert np.sum(a["positions"]) == 0
		break

	natoms = {'positions':[[0,0,0]]}
	sys.add_atoms(natoms)
	assert sys.natoms == 17

def test_embed():
	cu = bulk('Cu')
	sys = pc.System()
	sys.box = np.array(cu.cell)
	sys.atoms = {"positions": cu.positions}
	sys.embed_in_cubic_box()
	assert np.abs(sys.box[0][0] - 25.5265548) < 1E-5
