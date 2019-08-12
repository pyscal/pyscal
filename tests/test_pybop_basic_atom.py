import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs

def test_this_file():
    """
    Just a simple function to test if everything works
    """
    assert 1 == 1

def test_get_x():
    atom = pc.Atom(pos=[0,0,0], id=1)
    assert atom.get_x() == [0,0,0]
    atom.set_id(5)
    assert atom.get_id() == 5

def test_set_x():
    atom = pc.Atom(pos=[0,0,0], id=1)
    newx = [0,0,2]
    atom.set_x(newx)
    assert atom.get_x() == newx

def test_set_type():
    atom = pc.Atom(pos=[0,0,0], id=1, type=1)
    assert atom.get_type() == 1
    atom.set_type(2)
    assert atom.get_type() == 2

def test_set_solid():
    atom = pc.Atom(pos=[0,0,0], id=1, type=1)
    atom.set_solid(1)
    assert atom.get_solid() == 1
    with pytest.raises(ValueError):
        atom.set_solid(2)

def test_neighbors():
    atom = pc.Atom(pos=[0,0,0], id=1)
    atom.set_neighbors([1,2])
    assert atom.get_neighbors() == [1,2]
    assert atom.get_coordination() == 2
    atom.set_neighborweights([0.6, 0.4])
    assert atom.get_neighborweights() == [0.6, 0.4]