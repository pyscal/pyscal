"""
pyscal module for pickling of System and Atom objects. This module contains
pure python equivalents of C++/python hybrid classes of the same name. The
pure python classes can store the same information and can be pickled to save
the data.

Update : Atom objects are now pure - hence no pickling module is needed.
"""

import os
import numpy as np
import warnings
import pyscal.catom as pc

# Atom pickling support
class pickleAtom:
    """
    A pure python replica of Atom class for pickling information
    """
    def __init__(self):

        self.created = 1
        self.pos = [0,0,0]
        self.cluster = 0
        self.neighbors = []
        self.coordination = []
        self.neighbor_weights = []
        self.bonds = 0
        self.id = 0
        self.condition = 0
        self.solid = 0
        self.surface = 0
        self.largest_cluster = 0
        self.structure = 0
        self.loc = 0
        self.allq = []
        self.allaq = []
        self.type = 0
        self.custom = {}
        self.volume = 0
        self.avg_volume = 0
        self.face_vertices = []
        self.face_perimeters = []
        self.vertex_numbers = []
        self.vertex_vectors = []
        self.edge_lengths = []
        self.vorovector = []
        self.avg_connection = []

def pickle_atom(atom):
    """
    Get a python atom object and convert it to a picklable atom.
    Parameters
    ----------
    None
    Returns
    -------
    patom : picklable Atom object
    """
    patom = pickleAtom()

    patom.pos = atom.pos
    patom.cluster = atom.cluster
    patom.neighbors = atom.neighbors
    patom.coordination = atom.coordination
    patom.neighbor_weights = atom.neighbor_weights
    patom.bonds = atom.bonds
    patom.id = atom.id
    patom.condition = atom.condition
    patom.solid = atom.solid
    patom.surface = atom.surface
    patom.largest_cluster = atom.largest_cluster
    patom.structure = atom.structure
    patom.loc = atom.loc
    patom.allq = atom.allq
    patom.allaq = atom.allaq
    patom.type = atom.type
    patom.custom = atom.custom
    patom.volume = atom.volume
    patom.avg_volume = atom.avg_volume
    patom.face_vertices = atom.face_vertices
    patom.face_perimeters = atom.face_perimeters
    patom.vertex_numbers = atom.vertex_numbers
    patom.vertex_vectors = atom.vertex_vectors
    patom.edge_lengths = atom.edge_lengths
    patom.vorovector = atom.vorovector
    patom.avg_connection = atom.avg_connection

    return patom


def unpickle_atom(atom):
    """
    Get a picklable atom and convert it to a catom
    Parameters
    ----------
    patom : picklable Atom
    Returns
    -------
    catom : A c++ Atom object
    """
    patom = pc.Atom()

    patom.pos = atom.pos
    patom.cluster = atom.cluster
    patom.neighbors = atom.neighbors
    patom.coordination = atom.coordination
    patom.neighbor_weights = atom.neighbor_weights
    patom.bonds = atom.bonds
    patom.id = atom.id
    patom.condition = atom.condition
    patom.solid = atom.solid
    patom.surface = atom.surface
    patom.largest_cluster = atom.largest_cluster
    patom.structure = atom.structure
    patom.loc = atom.loc
    patom.allq = atom.allq
    patom.allaq = atom.allaq
    patom.type = atom.type
    patom.custom = atom.custom
    patom.volume = atom.volume
    patom.avg_volume = atom.avg_volume
    patom.face_vertices = atom.face_vertices
    patom.face_perimeters = atom.face_perimeters
    patom.vertex_numbers = atom.vertex_numbers
    patom.vertex_vectors = atom.vertex_vectors
    patom.edge_lengths = atom.edge_lengths
    patom.vorovector = atom.vorovector
    patom.avg_connection = atom.avg_connection

    return patom


def write_atoms(file, atoms):
    """
    Get a list of python atoms and write them to the given output file

    Parameters
    ----------
    file : string
        output file
    atoms : list of Atoms
        list of atom objects
    Returns
    -------
    None
    """
    patoms = np.array([pickle_atom(atom) for atom in atoms])
    np.save(file, patoms)

def read_atoms(file):
    """
    Read pickled atoms from a file and return catoms
    Parameters
    ----------
    file : string
        output file
    Returns
    -------
    atoms : list of atoms
        list of python atoms
    """
    patoms = np.load(file, allow_pickle=True).flatten()
    atoms = np.array([unpickle_atom(atom) for atom in patoms])
    return atoms


#System pickling
class pickleSystem:
    """
    A picklable system class
    """
    def __init__(self):
        self.created = True
        self.indicators = None
        self.atoms = None
        self.box = None
        self.rot = None
