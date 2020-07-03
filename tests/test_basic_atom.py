import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_this_file():
    """
    Just a simple function to test if everything works
    """
    assert 1 == 1

def test_get_x():
    atom = pc.Atom(pos=[0,0,0], id=1)
    assert atom.pos == [0,0,0]
    atom.id = 5
    assert atom.id == 5

def test_set_x():
    atom = pc.Atom(pos=[0,0,0], id=1)
    newx = [0,0,2]
    atom.pos = newx
    assert atom.pos == newx

def test_set_type():
    atom = pc.Atom(pos=[0,0,0], id=1, type=1)
    assert atom.type == 1
    atom.type = 2
    assert atom.type == 2

def test_set_solid():
    atom = pc.Atom(pos=[0,0,0], id=1, type=1)
    atom.solid = 1
    assert atom.solid == 1


def test_neighbors():
    atom = pc.Atom(pos=[0,0,0], id=1)
    atom.neighbors = [1,2]
    assert atom.neighbors == [1,2]
    assert atom.coordination == 2
    atom.neighbor_weights = [0.6, 0.4]
    assert atom.neighbor_weights == [0.6, 0.4]
