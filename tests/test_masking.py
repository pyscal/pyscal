import pyscal.core as pc
import os
import pyscal.crystal_structures as pcs
import numpy as np

def test_masking():
	atoms, box = pcs.make_crystal(structure='b2', repetitions=(6,6,6), lattice_constant=3.127)
	sys = pc.System()
	sys.box = box
	sys.atoms = atoms
	masks = [True if sys.types[x]==1 else False for x in range(sys.natoms)]
	masks2 = [False if sys.types[x]==1 else True for x in range(sys.natoms)]
	sys.apply_mask(masks2, mask_type="secondary")
	sys.apply_mask(masks, mask_type="primary")
	sys.find_neighbors(method="cutoff", cutoff=3.6)
	assert len(sys.neighbors[0]) == 0
	assert len(sys.neighbors[1]) == 7

	sys.remove_mask(mask_type="all")
	sys.find_neighbors(method="cutoff", cutoff=3.6)
	assert len(sys.neighbors[0]) == 14
	assert len(sys.neighbors[1]) == 14
