import pyscal.core as pc
import os
import pyscal.crystal_structures as pcs
import numpy as np
from ase.build import bulk

def test_system_init():
	sys = pc.System()
	sys.read_inputfile("tests/files/conf.dump", customkeys=["vx", "vy"])
	assert sys.natoms == 500
	assert sys.atoms.vx[0] == '0.0394436'

	sys.read_inputfile("tests/files/conf.dump.gz")
	assert sys.natoms == 500

	sys.read_inputfile("tests/files/conf.bcc.scaled.dump")
	assert sys.natoms == 2

	sys.read_inputfile("tests/files/POSCAR", format="poscar")
	assert sys.natoms == 42

	sys = pc.System()
	sys.read_inputfile("tests/files/conf.dump", customkeys=["vx", "vy"])
	sys.to_file("test.dump", customkeys=["vx"])

	sys = pc.System()
	sys.read_inputfile("test.dump", customkeys=["vx"])
	assert sys.natoms == 500
	assert sys.atoms.vx[0] == '0.0394436'

	sys = pc.System()
	sys.read_inputfile("tests/files/conf.dump", customkeys=["vx", "vy"])
	aseobj = sys.to_ase(species=["Au"])
	assert np.sum(sys.atoms.positions[0]-aseobj.positions[0]) < 1E-5