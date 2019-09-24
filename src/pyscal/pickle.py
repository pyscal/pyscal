"""
pyscal module for pickling of System and Atom objects. This module contains
pure python equivalents of C++/python hybrid classes of the same name. The
pure python classes can store the same information and can be pickled to save
the data.

Update : Atom objects are now pure - hence no pickling module is needed.
"""

import pyscal.ccore as pc
#import pyscal.core as pcc
import pyscal.traj_process as ptp
import os
import numpy as np
import warnings

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


#def write_systems(file, systems):
    """
    #Write an array of systems to file

    #Parameters
    #----------
    #file : string
    #    name of output file

    #systems : list of system objects

    #Returns
    #-------
    #None

    """
    #psystems = [sys.prepare_pickle() for sys in systems]
    #np.save(file, psystems)

#def read_systems(file):
    """
    #Unpickle systems from file

    #Parameters
    #----------
    #file : string
    #    name of input file

    #Returns
    #-------
    #None
    """

    #psystems = np.load(file, allow_pickle=True).flatten()
    #systems = np.array([unpickle_system(system) for system in psystems])
    #return systems
