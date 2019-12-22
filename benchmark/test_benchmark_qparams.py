import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

sys = pc.System()
atoms, box = pcs.make_crystal('fcc', lattice_constant=4, repetitions=[10,10,10])
sys.atoms = atoms
sys.box = box
sys.find_neighbors(method='cutoff', cutoff=0)

def test_q2(benchmark):
    benchmark(sys.calculate_q, 2)

def test_q3(benchmark):
    benchmark(sys.calculate_q, 3)

def test_q4(benchmark):
    benchmark(sys.calculate_q, 4)

def test_q5(benchmark):
    benchmark(sys.calculate_q, 5)

def test_q6(benchmark):
    benchmark(sys.calculate_q, 6)

def test_q7(benchmark):
    benchmark(sys.calculate_q, 7)

def test_q8(benchmark):
    benchmark(sys.calculate_q, 8)

def test_q9(benchmark):
    benchmark(sys.calculate_q, 9)

def test_q10(benchmark):
    benchmark(sys.calculate_q, 10)

def test_q11(benchmark):
    benchmark(sys.calculate_q, 11)

def test_q12(benchmark):
    benchmark(sys.calculate_q, 12)
