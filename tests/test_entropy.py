import pytest
import os
import numpy as np
import pyscal.core as pc


def test_entropy():
	sys = pc.System()
	sys.read_inputfile("tests/conf.fcc.Al.dump")
	sys.find_neighbors(method="cutoff", cutoff=0)
	lat = (sys.box[0][1]-sys.box[0][0])/5
	sys.calculate_entropy(1.4*lat, ra=0.9*lat, averaged=True)
	atoms = sys.atoms
	solid_entropy = [atom.entropy for atom in atoms]
	solid_avg_entropy = [atom.avg_entropy for atom in atoms]
	assert np.abs(np.mean(solid_entropy) + 3.47249) < 0.001
	assert np.abs(np.mean(solid_avg_entropy) + 3.47254) < 0.001