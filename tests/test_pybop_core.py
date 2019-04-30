import pytest
import os,sys,inspect

import pybop.core as pc

atom = pc.Atom(pos=[0,0,0], id=1)
natom = pc.Atom(pos=[0,0,1], id=2)
sys = pc.System()
sys.read_inputfile("conf.dump")

def test_this_file():
    """
    Just a simple function to test if everything works
    """
    assert 1 == 1

def test_get_x():
    assert atom.get_x() == [0,0,0]
    atom.set_id(5)
    assert atom.get_id() == 5


def test_set_x():
    newx = [0,0,2]
    atom.set_x(newx)
    assert atom.get_x() == newx

def test_neighbors():
    atom.set_neighbors([1,2])
    assert atom.get_neighbors() == [1,2]
    assert atom.get_coordination() == 2
    atom.set_neighborweights([0.6, 0.4])
    assert atom.get_neighborweights() == [0.6, 0.4]

def test_basic_system():
    #basic system tests
    sys = pc.System()
    sys.set_box([[0,1],[0,1],[0,1]])
    assert sys.get_box() == [[0,1],[0,1],[0,1]]
    #sys.read_inputfile("conf.dump")


def test_neighbors_system():
    #first test the cutoff method
    pass


