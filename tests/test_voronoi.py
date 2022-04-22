import pyscal.core as pc
import os
import pyscal.crystal_structures as pcs
import numpy as np
from ase.build import bulk

def test_voronoi_props():
	nx = 5
	satoms, box = pcs.make_crystal(structure='bcc', 
    	lattice_constant=3.127, 
    	repetitions=(nx,nx,nx),)
	sys = pc.System()
	sys.box = box
	sys.atoms = satoms
	sys.find_neighbors(method="voronoi")

	assert sys.atom.voronoi.vertex.numbers[0][0] == 6
	assert (sys.atom.voronoi.vertex.positions[0][0][0]+1.5635 < 1E-4)
	assert (sys.atom.voronoi.volume[0]-15.288104691499992 < 1E-5)
	assert (sys.atom.voronoi.face.perimeters[0][0]-6.6333687143110005 < 1E-5)
	assert sys.atom.voronoi.face.vertices[0][0] == 6
