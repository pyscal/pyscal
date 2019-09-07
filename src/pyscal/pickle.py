import pyscal.ccore as pc
import pyscal.core as pcc
import pyscal.traj_process as ptp
import os
import numpy as np
import warnings

"""
Pickling module for pyscal
"""

class Atom:
    """
    A pure python replica of Atom class for pickling information
    """
    def __init__(self):
        self.created = 1

        self.pos = None
        self.solid = None
        self.structure = None
        self.cluster = None
        self.neighbors = None
        self.neighborweights = None
        self.allq = None
        self.allaq = None
        self.id = None
        self.loc = None
        self.type = None
        self.vorovector = None
        self.volume = None
        self.avgvolume = None
        self.facevertices = None

def get_picklable_atom(atom):
    """
    Get a python atom object and convert it to a picklable atom.

    Parameters
    ----------
    None

    Returns
    -------
    patom : picklable Atom object

    """
    patom = Atom()

    patom.pos = atom.get_x()
    patom.solid = atom.get_solid()
    patom.structure = atom.get_structure()
    patom.cluster = atom.get_cluster()
    patom.neighbors = atom.get_neighbors()
    patom.neighborweights = atom.get_neighborweights()
    patom.allq = atom.get_allq()
    patom.allaq = atom.get_allaq()
    patom.id = atom.get_id()
    patom.loc = atom.get_loc()
    patom.type = atom.get_type()
    patom.vorovector = atom.get_vorovector()
    patom.volume = atom.get_volume()
    patom.avgvolume = atom.get_avgvolume()
    patom.facevertices = atom.get_facevertices()

    return patom


def convert_picklable_atom(patom):
    """
    Get a picklable atom and convert it to a catom

    Parameters
    ----------
    patom : picklable Atom

    Returns
    -------
    catom : A c++ Atom object
    """
    atomc = pcc.Atom()
    atomc.set_pos(patom.pos)
    atomc.set_solid(patom.solid)
    atomc.set_structure(patom.structure)
    atomc.set_cluster(patom.cluster)
    atomc.set_neighbors(patom.neighbors)
    atomc.set_neighborweights(patom.neighborweights)
    atomc.set_allq(patom.allq)
    atomc.set_allaq(patom.allaq)
    atomc.set_id(patom.id)
    atomc.set_loc(patom.loc)
    atomc.set_type(patom.type)
    atomc.set_vorovector(patom.vorovector)
    atomc.set_volume(patom.volume)
    atomc.set_avgvolume(patom.avgvolume)
    atomc.set_facevertices(patom.facevertices)
    return atomc    


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
    patoms = np.array([get_picklable_atom(atom) for atom in atoms])
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
    patoms = np.load(file, allow_pickle=True)
    atoms = np.array([convert_picklable_atom(atom) for atom in patoms])
    return atoms

#System pickling
class System:
    """
    A picklable system class
    """
    def __init__(self):
        self.created = True

    self.atoms = None
    self.boxdims = None
    self.alpha = None
    self.rot = None
    

