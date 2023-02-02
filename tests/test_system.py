import pyscal.core as pc
import os
import numpy as np
from ase.build import bulk
from pyscal.atoms import Atoms
from pyscal.crystal_structures import Structure


def test_system_init():
	sys = Structure().lattice.bcc(repetitions = [4, 4, 4], lattice_constant=3.127)
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
	sys = Structure().lattice.bcc(repetitions = [2, 2, 2], lattice_constant=3.127)

	assert sys.natoms == 16
	assert len(sys.atoms['positions']) == 432

	for a in sys.iter_atoms():
		assert np.sum(a["positions"]) == 0
		break

	natoms = {'positions':[[0,0,0]]}
	sys.atoms += natoms
	assert sys.natoms == 17

def test_embed():
	cu = bulk('Cu')
	sys = pc.System()
	sys.box = np.array(cu.cell)
	sys.atoms = Atoms({"positions": cu.positions})
	sys.embed_in_cubic_box()
	assert np.abs(sys.box[0][0] - 15.315932880500618) < 1E-5

def test_distance():
	sys = Structure().lattice.bcc(repetitions = [2, 2, 2], lattice_constant=3.127)
	dist = sys.get_distance([0.0, 0.0, 0.0], [1.5635, 1.5635, 1.5635])
	assert np.abs(dist - 2.708061437633939) < 1E-5	

def test_composition():
	sys = Structure().lattice.l12(repetitions = [2, 2, 2], lattice_constant=3.127)
	c = sys.concentration
	assert c['1'] == 8
	assert c['2'] == 24

	c = sys.composition
	assert c['1'] == 8
	assert c['2'] == 24

def test_volume():
	sys = Structure().lattice.fcc(repetitions = [10, 10, 10])
	assert sys.volume == 1000