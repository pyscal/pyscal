import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

sys = pc.System()
atoms, box = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[4,4,4])
sys.atoms = atoms
sys.box = box
sys.find_neighbors(method='cutoff', cutoff=0)

def test_chiparams(benchmark):
    benchmark(sys.calculate_chiparams)

def test_disorder(benchmark):
    sys.calculate_q(6)
    benchmark(sys.calculate_disorder)

def test_find_solids(benchmark):
    benchmark(sys.find_solids, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)
