"""
This file will run quick benchmark tests for all the major functions
"""
import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

sys = pc.System()
atoms, box = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[13,13,13])
sys.atoms = atoms
sys.box = box

def test_benchmark_neigh_cutoff(benchmark):
    #firstly benchmark neighbor methods
    benchmark(sys.find_neighbors, method='cutoff', cutoff=3.7)

def test_benchmark_neigh_adapt(benchmark):
    benchmark(sys.find_neighbors, method='cutoff', cutoff=0)

def test_benchmark_neigh_sann(benchmark):
    benchmark(sys.find_neighbors, method='cutoff', cutoff='sann')

def test_benchmark_neigh_voro(benchmark):
    benchmark(sys.find_neighbors, method='voronoi')
