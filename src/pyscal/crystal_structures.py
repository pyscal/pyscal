"""
pyscal module for creating crystal structures.
"""

import numpy as np
import warnings
import os
from functools import partial
from functools import update_wrapper


from pyscal.attributes import read_yaml
from pyscal.core import System
from pyscal.structure_creator import make_crystal

structures = read_yaml(os.path.join(os.path.dirname(__file__), "data/structure_data.yaml"))
elements = read_yaml(os.path.join(os.path.dirname(__file__), "data/element_data.yaml"))


#wrapper methods
def structure_creator(structure, 
             lattice_constant = 1.00, 
             repetitions = None, 
             ca_ratio = 1.633, 
             noise = 0,
             element=None):
    """
    Create a crystal structure and return it as a System object.

    Parameters
    ----------
    structure : {'sc', 'bcc', 'fcc', 'hcp', 'diamond', 'a15' or 'l12'}
        type of the crystal structure

    lattice_constant : float, optional
        lattice constant of the crystal structure, default 1

    repetitions : list of ints of len 3, optional
        of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
        default `[1, 1, 1]`.

    ca_ratio : float, optional
        ratio of c/a for hcp structures, default 1.633

    noise : float, optional
        If provided add normally distributed noise with standard deviation `noise` to the atomic positions.
    
    element : string, optional
        The chemical element

    Returns
    -------
    System: pyscal System
        system will be populated with given atoms and simulation box

    Examples
    --------
    >>> sys = structure_creator('bcc', lattice_constant=3.48, repetitions=[2,2,2])
    """    
    return System.from_structure(structure, lattice_constant=lattice_constant,
             repetitions=repetitions, ca_ratio=ca_ratio,
             noise=noise, element=element)

class Structure:
    """
    A class for structure creation
    
    Attributes
    ----------
    element: 
        create elementary structures
    lattice:
        create structures by specifying lattice
    """
    def __init__(self):
        #create by element name
        self.element = ElementCreator(elements)
        #create by lattice name
        self.lattice = LatticeCreator(structures)
        #create a general structure
        self.custom = general_lattice
        #complete dict
        self._structure_dict = structures

    def structure_dict(self, structure):
        return self._structure_dict[structure]

class ElementCreator:
    """
    Create an elementary structure
    """
    def __init__(self, element_dict):
        self._element_dict = element_dict
    
    def __dir__(self):
        return list(self._element_dict.keys())
    
    def __getattr__(self, key):
        #this is the element based creater
        if key in self._element_dict.keys():
            structure = self._element_dict[key]['structure']
            pfunc = partial(structure_creator, structure,
                        lattice_constant=self._element_dict[key]['lattice_constant'],
                        element = key)
            update_wrapper(pfunc, structure_creator)
            return pfunc

class LatticeCreator(ElementCreator):
    """
    Create a lattice
    """
    def __getattr__(self, key):
        if key in self._element_dict.keys():
            pfunc = partial(structure_creator, key)
            update_wrapper(pfunc, structure_creator)
            return pfunc


#general structure creator
def general_lattice(species, positions, 
    scaling_factors=[1.0, 1.0, 1.0],
    lattice_constant = 1.00, 
    repetitions = None, 
    noise = 0,
    element=None):
    """
    Create a general lattice structure.

    species: list
        list of per-atom species

    positions:
        list of relative positions positions of reach atom (between 0-1)

    scaling_fractors:
        factors with which the unit cell should be scaled, for example hcp could
        have [1,1.73, 1.63]. Default [1,1,1]

    lattice_constant : float, optional
        lattice constant of the crystal structure, default 1

    repetitions : list of ints of len 3, optional
        of type `[nx, ny, nz]`, repetions of the unit cell in x, y and z directions.
        default `[1, 1, 1]`.

    noise : float, optional
        If provided add normally distributed noise with standard deviation `noise` to the atomic positions.
    
    element : string, optional
        The chemical element
    """
    if not (len(species) == len(positions)):
        raise ValueError("Species and positions should have same length!")

    sdict = {"custom":
                {"natoms": len(positions),
                 "species": species,
                 "scaling_factors": scaling_factors,
                 "positions": positions}
            }

    atoms, box = make_crystal("custom", lattice_constant=lattice_constant,
        repetitions=repetitions, noise=noise, element=element,
        structures=sdict)

    sys = System()
    sys.box = box
    sys.atoms = atoms
    sys.atoms._lattice = None
    sys.atoms._lattice_constant = lattice_constant
    sys._structure_dict = sdict["custom"]
    return sys

        

    
    

