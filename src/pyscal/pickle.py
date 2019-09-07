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


def unpickle_atom(patom):
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
class System:
    """
    A picklable system class
    """
    def __init__(self):
        self.created = True

        self.indicators = None

        self.atoms = None
        self.boxdims = None
        self.rot = None


def unpickle_system(psys):
    """
    Unpickle a picklable system to proper system
    
    Parameters
    ----------
    psys : pickled system

    Returns
    -------
    sys : unpickled system

    """
    #create a system instance
    sys = pcc.System()

    #set indicators
    sys.set_indicators(psys.indicators)

    #set atoms
    hatoms = [unpickle_atom(atom) for atom in psys.atoms]
    #convert to catoms
    catoms = [sys.copy_atom_to_catom(atom) for atom in hatoms]

    #get box
    boxdims = psys.boxdims

    #if triclinic, get those
    if psys.indicators[6] == 1:
        rot = psys.rot
        rotinv = np.linalg.inv(rot)
        sys.assign_triclinic_params(rot, rotinv)

    #assign atoms and box
    sys.reassign_particles(catoms, boxdims)

    return sys


def write_systems(file, systems):
    """
    Write an array of systems to file

    Parameters
    ----------
    file : string
        name of output file

    systems : list of system objects

    Returns
    -------
    None

    """
    psystems = [sys.prepare_pickle() for sys in systems]
    np.save(file, psystems)

def read_systems(file):
    """
    Unpickle systems from file

    Parameters
    ----------
    file : string
        name of input file

    Returns
    -------
    None
    """

    psystems = np.load(file, allow_pickle=True).flatten()
    systems = np.array([unpickle_system(system) for system in psystems])
    return systems








    

