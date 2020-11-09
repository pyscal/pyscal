import numpy as np
import gzip
import pyscal.catom as pca
from ase import Atom, Atoms
import gzip
import io
import os
from ase.io import write, read
import pyscal.trajectory.ase as ptase

def read_snap(infile, compressed = False):
    """
    Function to read a POSCAR format.

    Parameters
    ----------
    infile : string
        name of the input file

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary, Default False

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    box : list of list of floats
        list of the type `[[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]]` where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = read_poscar('POSCAR')
    >>> atoms, box = read_poscar('POSCAR.gz')
    >>> atoms, box = read_poscar('POSCAR.dat', compressed=True)

    """
    aseobj = read(infile, format="vasp")
    atoms, box = ptase.read_snap(aseobj)
    return atoms, box


def write_snap(sys, outfile, comments="pyscal", species=None):
    """
    Function to read a POSCAR format.

    Parameters
    ----------
    outfile : string
        name of the input file


    """
    aseobj = ptase.convert_snap(sys, species=species)
    write(outfile, aseobj, format="vasp")


def split_snaps(**kwargs):
    raise NotImplementedError("split method for mdtraj is not implemented")

def convert_snap(**kwargs):
    raise NotImplementedError("convert method for mdtraj is not implemented")
