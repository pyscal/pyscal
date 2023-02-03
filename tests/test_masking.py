import pyscal.core as pc
import os
import pyscal.crystal_structures as pcs
import numpy as np

def test_masking():
	pass
	#atoms, box = pcs.make_crystal(structure='b2', repetitions=(6,6,6), lattice_constant=3.127)
	#sys = pc.System()
	#sys.box = box
	#sys.atoms = atoms

	#def _condition1(atom):
	#	if atom.species[0] == 1:
	#		return True
	#	return False

	#def _condition2(atom):
	#	if atom.species[0] == 1:
	#		return False
	#	return True

	#sys.apply_mask(condition=_condition1, mask_type="secondary")
	#sys.apply_mask(condition=_condition2, mask_type="primary")	
	#sys.find_neighbors(method="cutoff", cutoff=3.6)
	#assert len(sys.atoms.neighbors.index[0]) == 0
	#sys.remove_mask(mask_type="all")

	#sys.apply_mask(condition=_condition2, mask_type="primary")
	#sys.find_neighbors(method="cutoff", cutoff=3.6)
	#assert len(sys.atoms.neighbors.index[1]) == 7
	#sys.remove_mask(mask_type="all")

	#sys.find_neighbors(method="cutoff", cutoff=3.6)
	#assert len(sys.atoms.neighbors.index[0]) == 14
	#assert len(sys.atoms.neighbors.index[1]) == 14
