import pyscal.core as pc
import os
from pyscal.crystal_structures import Structure
import numpy as np
from ase.build import bulk

def test_voronoi_props():
    nx = 5
    sys = Structure().lattice.bcc(repetitions = [nx, nx, nx], lattice_constant=3.127)
    sys.find_neighbors(method="voronoi")

    assert sys.atom.voronoi.vertex.numbers[0][0] == 6
    assert (sys.atom.voronoi.vertex.positions[0][0][0]+1.5635 < 1E-4)
    assert (sys.atom.voronoi.volume[0]-15.288104691499992 < 1E-5)
    assert (sys.atom.voronoi.face.perimeters[0][0]-6.6333687143110005 < 1E-5)
    assert sys.atom.voronoi.face.vertices[0][0] == 6

def test_voronoi_vertices():
    nx = np.random.randint(1, 10)
    nverts = (nx**3*12)
    sys = Structure().lattice.bcc(repetitions = [nx, nx, nx])
    sys.find_neighbors(method='voronoi', cutoff=0.1)
    assert len(sys.atom.voronoi.vertex.unique_positions) == nverts